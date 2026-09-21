using Chromatography;                        // SeparationType (NOT CommonParameters.SeparationType, which is a string)
using Chromatography.RetentionTimePrediction;
using Omics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.FdrAnalysis
{
    /// <summary>
    /// An <see cref="IRetentionTimePredictor"/> that answers from a table built up-front in batches,
    /// so PEP never calls the model one peptide at a time.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Why this exists. Chronologer runs its TorchSharp model inside a process-wide lock, and PEP asked it
    /// for one peptide at a time from a 32-thread <c>Parallel.ForEach</c>. On a single low-resolution file
    /// that is ~99,000 locked calls, and a wall-clock profile put 66.99% of the whole task inside
    /// <c>ChronologerRetentionTimePredictor.PredictCore</c> with 64.36% of it blocked in
    /// <c>Monitor.Enter_Slowpath</c>. The model itself was ~5%; the queue to reach it was ~62%. Every core
    /// was busy and almost none of it was inference.
    /// </para>
    /// <para>
    /// mzLib already ships the cure: <see cref="IRetentionTimePredictor.PredictRetentionTimeEquivalents"/>
    /// takes the lock once per call and runs 2048-row forward passes. Warming through it turns those ~99,000
    /// lock acquisitions into ~43.
    /// </para>
    /// <para>
    /// This is a decorator rather than a subclass because
    /// <c>RetentionTimePredictor.PredictRetentionTimeEquivalent</c> is not virtual, so a derived class cannot
    /// intercept the single-peptide path. It does NOT own the predictor it wraps, and deliberately does not
    /// dispose it: MetaMorpheus shares one Chronologer through a process-lifetime <c>static Lazy</c>.
    /// </para>
    /// <para>
    /// A miss falls through to the wrapped predictor, so behaviour is unchanged for anything not warmed.
    /// </para>
    /// </remarks>
    internal sealed class PrewarmedRetentionTimePredictor : IRetentionTimePredictor
    {
        /// <summary>
        /// Peptides handed to the batched API at once. Its first phase allocates a 52-long <c>long[]</c> per
        /// peptide before any inference starts, so the whole corpus in one call would be gigabytes of encode
        /// buffer for no benefit.
        /// </summary>
        private const int WarmChunkSize = 50_000;

        private readonly IRetentionTimePredictor _inner;

        /// <summary>
        /// Keyed on full sequence, because that is what the prediction is a function of.
        /// </summary>
        /// <remarks>
        /// The value is the PAIR, not the value alone. <c>(value, reason)</c> is not "value XOR reason":
        /// a predictor in <c>UsePrimarySequence</c> mode reports a failure reason AND returns a usable
        /// sequence, so a successful prediction can legitimately arrive with a reason attached. Storing only
        /// the double would silently discard the diagnostic that <see cref="PsmData.HasHydrophobicity"/>
        /// exists to carry.
        /// </remarks>
        private readonly Dictionary<string, (double? Value, RetentionTimeFailureReason? Reason)> _byFullSequence;

        private PrewarmedRetentionTimePredictor(
            IRetentionTimePredictor inner,
            Dictionary<string, (double? Value, RetentionTimeFailureReason? Reason)> byFullSequence)
        {
            _inner = inner;
            _byFullSequence = byFullSequence;
        }

        /// <summary>Distinct full sequences warmed. Logged so the batching ratio stays visible.</summary>
        internal int WarmedSequenceCount => _byFullSequence.Count;

        /// <summary>
        /// Predicts every distinct peptidoform in <paramref name="peptides"/> in batches and returns a
        /// predictor that serves those answers without touching the model again.
        /// </summary>
        /// <remarks>
        /// Takes peptides rather than PSMs on purpose: what this type caches is a function of the
        /// peptidoform alone, and keeping spectral-match types out of it makes the caching behaviour
        /// testable without constructing a search result.
        /// </remarks>
        internal static PrewarmedRetentionTimePredictor Warm(
            IRetentionTimePredictor inner, IEnumerable<IRetentionPredictable> peptides, int maxThreads)
        {
            if (inner == null)
            {
                throw new ArgumentNullException(nameof(inner));
            }

            // One representative per distinct full sequence. Chronologer's batched override does not
            // deduplicate -- it encodes every repeated row -- so collapsing here is the caller's job.
            Dictionary<string, IRetentionPredictable> representatives = new();
            foreach (IRetentionPredictable peptide in peptides ?? Enumerable.Empty<IRetentionPredictable>())
            {
                string fullSequence = peptide?.FullSequence;
                if (!string.IsNullOrEmpty(fullSequence))
                {
                    representatives.TryAdd(fullSequence, peptide);
                }
            }

            Dictionary<string, (double?, RetentionTimeFailureReason?)> warmed = new(representatives.Count);
            foreach (IRetentionPredictable[] chunk in Chunk(representatives.Values, WarmChunkSize))
            {
                // Two things about this loop. maxThreads parallelises the batched call's CPU-side
                // format-and-encode phase only -- inference is serialised under the model lock either way.
                // And results are keyed back by the tuple's OWN peptide rather than by position, because
                // IRetentionTimePredictor does not promise the result order matches the input order.
                foreach ((double? value, IRetentionPredictable peptide, RetentionTimeFailureReason? reason)
                         in inner.PredictRetentionTimeEquivalents(chunk, maxThreads))
                {
                    // PredictionError is the bucket a blanket catch puts transient native faults into
                    // (mzLib #1075 fixed one that only appeared in high-volume runs). Every other reason is
                    // a pure function of the peptidoform and is safe to reuse; this one is not, so leave it
                    // out and let the miss path retry it against the model.
                    if (reason == RetentionTimeFailureReason.PredictionError)
                    {
                        continue;
                    }

                    string fullSequence = peptide?.FullSequence;
                    if (!string.IsNullOrEmpty(fullSequence))
                    {
                        warmed[fullSequence] = (value, reason);
                    }
                }
            }

            return new PrewarmedRetentionTimePredictor(inner, warmed);
        }

        /// <summary>Splits the warm set into batches small enough that the encode phase stays bounded.</summary>
        private static IEnumerable<IRetentionPredictable[]> Chunk(IEnumerable<IRetentionPredictable> source, int size)
        {
            List<IRetentionPredictable> buffer = new(size);
            foreach (IRetentionPredictable item in source)
            {
                buffer.Add(item);
                if (buffer.Count == size)
                {
                    yield return buffer.ToArray();
                    buffer.Clear();
                }
            }

            if (buffer.Count > 0)
            {
                yield return buffer.ToArray();
            }
        }

        public string PredictorName => _inner.PredictorName;

        public SeparationType SeparationType => _inner.SeparationType;

        public double? PredictRetentionTimeEquivalent(
            IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
        {
            // FullSequence is annotated non-nullable but the interface's own documentation says it may be
            // absent, and implementations defend against that, so do not assume it here either.
            string fullSequence = peptide?.FullSequence;
            if (!string.IsNullOrEmpty(fullSequence)
                && _byFullSequence.TryGetValue(fullSequence, out (double? Value, RetentionTimeFailureReason? Reason) warmed))
            {
                failureReason = warmed.Reason;
                return warmed.Value;
            }

            return _inner.PredictRetentionTimeEquivalent(peptide, out failureReason);
        }

        public IReadOnlyList<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)>
            PredictRetentionTimeEquivalents(IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1)
        {
            List<(double?, IRetentionPredictable, RetentionTimeFailureReason?)> results = new();
            List<IRetentionPredictable> misses = null;

            foreach (IRetentionPredictable peptide in peptides ?? Enumerable.Empty<IRetentionPredictable>())
            {
                string fullSequence = peptide?.FullSequence;
                if (!string.IsNullOrEmpty(fullSequence)
                    && _byFullSequence.TryGetValue(fullSequence, out (double? Value, RetentionTimeFailureReason? Reason) warmed))
                {
                    results.Add((warmed.Value, peptide, warmed.Reason));
                }
                else
                {
                    (misses ??= new List<IRetentionPredictable>()).Add(peptide);
                }
            }

            if (misses != null)
            {
                results.AddRange(_inner.PredictRetentionTimeEquivalents(misses, maxThreads));
            }

            return results;
        }

        public string GetFormattedSequence(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
            => _inner.GetFormattedSequence(peptide, out failureReason);

        /// <summary>
        /// Deliberately does nothing. This type borrows its inner predictor and must not dispose it -- the
        /// Chronologer instance is shared process-wide, and its own Dispose does not take the model lock,
        /// so disposing it while another thread is predicting is a race.
        /// </summary>
        public void Dispose()
        {
        }
    }
}
