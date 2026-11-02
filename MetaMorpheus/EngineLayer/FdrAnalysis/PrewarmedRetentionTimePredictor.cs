using Chromatography;                        // SeparationType (NOT CommonParameters.SeparationType, which is a string)
using Chromatography.RetentionTimePrediction;
using Omics;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading;

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
    /// A miss falls through to the wrapped predictor, so behaviour is unchanged for anything not warmed. Each
    /// missed peptidoform is asked about once: the answer is memoized, and the calls are serialised because
    /// Koina's <c>RetentionTimeModel</c> documents its instances as not thread-safe while the misses come from
    /// PEP's parallel loop.
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
        /// Stores whatever pair the batched call returned. Do not rely on a warmed reason matching what the
        /// single-peptide path would have said: Chronologer's batched override (mzLib 1.0.591) returns a null
        /// reason for every successful row, so a <c>UsePrimarySequence</c> reason the single path keeps is
        /// already gone. Nothing reads the reason today -- both consumers pass <c>out _</c>, and
        /// <see cref="PsmData.HasHydrophobicity"/> comes from whether a value exists.
        /// </remarks>
        private readonly Dictionary<string, (double? Value, RetentionTimeFailureReason? Reason)> _byFullSequence;

        /// <summary>
        /// Answers for peptidoforms that were not warmed, recorded the first time the wrapped predictor is
        /// asked, so a miss costs one call to the model rather than one per hypothesis.
        /// </summary>
        private readonly ConcurrentDictionary<string, (double? Value, RetentionTimeFailureReason? Reason)> _missed = new();

        /// <summary>Serialises every call that reaches the wrapped predictor after warming.</summary>
        private readonly object _innerLock = new();

        private int _missCount;

        private PrewarmedRetentionTimePredictor(
            IRetentionTimePredictor inner,
            Dictionary<string, (double? Value, RetentionTimeFailureReason? Reason)> byFullSequence,
            int failedChunkCount,
            string firstChunkFailure)
        {
            _inner = inner;
            _byFullSequence = byFullSequence;
            FailedChunkCount = failedChunkCount;
            FirstChunkFailure = firstChunkFailure;
        }

        /// <summary>Distinct full sequences warmed. Logged so the batching ratio stays visible.</summary>
        internal int WarmedSequenceCount => _byFullSequence.Count;

        /// <summary>
        /// Calls that reached the wrapped predictor because the answer was not warmed. Logged beside
        /// <see cref="WarmedSequenceCount"/>, so a regression back to one-at-a-time prediction shows up in the log.
        /// </summary>
        internal int MissCount => Volatile.Read(ref _missCount);

        /// <summary>Batched calls that threw during warming. Their peptidoforms fall through to the miss path.</summary>
        internal int FailedChunkCount { get; }

        /// <summary>The message of the first batched call that threw, or null when none did.</summary>
        internal string FirstChunkFailure { get; }

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
            int failedChunkCount = 0;
            string firstChunkFailure = null;
            foreach (IRetentionPredictable[] chunk in Chunk(representatives.Values, WarmChunkSize))
            {
                // Two things about this loop. maxThreads parallelises the batched call's CPU-side
                // format-and-encode phase only -- inference is serialised under the model lock either way.
                // And results are keyed back by the tuple's OWN peptide rather than by position, because
                // IRetentionTimePredictor does not promise the result order matches the input order.
                IReadOnlyList<(double?, IRetentionPredictable, RetentionTimeFailureReason?)> batch;
                try
                {
                    batch = inner.PredictRetentionTimeEquivalents(chunk, maxThreads);
                }
                catch (Exception e)
                {
                    // The single-peptide path turns a thrown exception into one PredictionError; Chronologer's
                    // batched override has no such catch, so without this one fault would fail the whole task.
                    // Leave the chunk unwarmed and let its peptidoforms take the miss path instead.
                    failedChunkCount++;
                    firstChunkFailure ??= e.GetBaseException().Message;
                    continue;
                }

                foreach ((double? value, IRetentionPredictable peptide, RetentionTimeFailureReason? reason) in batch)
                {
                    // Koina marks every row of an HTTP batch that failed as PredictionError, which may well
                    // succeed on a retry, so leave those out and give each one a single, memoized retry on the
                    // miss path. Chronologer returns PredictionError only for format and encode failures, which
                    // its single path rejects before taking the model lock, so the retry is cheap there.
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

            return new PrewarmedRetentionTimePredictor(inner, warmed, failedChunkCount, firstChunkFailure);
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

            if (string.IsNullOrEmpty(fullSequence))
            {
                lock (_innerLock)
                {
                    Interlocked.Increment(ref _missCount);
                    return _inner.PredictRetentionTimeEquivalent(peptide, out failureReason);
                }
            }

            (double? Value, RetentionTimeFailureReason? Reason) missed = PredictMiss(peptide, fullSequence);
            failureReason = missed.Reason;
            return missed.Value;
        }

        /// <summary>
        /// Asks the wrapped predictor about a peptidoform that was not warmed, once, and remembers the answer.
        /// </summary>
        private (double? Value, RetentionTimeFailureReason? Reason) PredictMiss(IRetentionPredictable peptide, string fullSequence)
        {
            if (_missed.TryGetValue(fullSequence, out var known))
            {
                return known;
            }

            lock (_innerLock)
            {
                if (_missed.TryGetValue(fullSequence, out known))
                {
                    return known;
                }

                Interlocked.Increment(ref _missCount);
                double? value = _inner.PredictRetentionTimeEquivalent(peptide, out RetentionTimeFailureReason? reason);
                _missed[fullSequence] = (value, reason);
                return (value, reason);
            }
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
                else if (!string.IsNullOrEmpty(fullSequence) && _missed.TryGetValue(fullSequence, out var known))
                {
                    results.Add((known.Value, peptide, known.Reason));
                }
                else
                {
                    (misses ??= new List<IRetentionPredictable>()).Add(peptide);
                }
            }

            if (misses != null)
            {
                lock (_innerLock)
                {
                    Interlocked.Add(ref _missCount, misses.Count);
                    foreach ((double? value, IRetentionPredictable peptide, RetentionTimeFailureReason? reason)
                             in _inner.PredictRetentionTimeEquivalents(misses, maxThreads))
                    {
                        string fullSequence = peptide?.FullSequence;
                        if (!string.IsNullOrEmpty(fullSequence))
                        {
                            _missed.TryAdd(fullSequence, (value, reason));
                        }
                        results.Add((value, peptide, reason));
                    }
                }
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
