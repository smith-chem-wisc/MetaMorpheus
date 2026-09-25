using Chromatography;
using Chromatography.RetentionTimePrediction;
using EngineLayer.FdrAnalysis;
using NUnit.Framework;
using Omics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal class PrewarmedRetentionTimePredictorTest
    {
        /// <summary>
        /// Minimal <see cref="IRetentionPredictable"/>; the predictor only ever keys on FullSequence.
        /// </summary>
        private sealed class StubPeptide : IRetentionPredictable
        {
            public string BaseSequence { get; init; } = string.Empty;
            public string FullSequence { get; init; } = string.Empty;
            public string FullSequenceWithMassShifts { get; init; } = string.Empty;
            public double MonoisotopicMass { get; init; } = 1000;
        }

        /// <summary>
        /// Returns whatever it is told to, and counts how many times each entry point was used, so a test
        /// can assert the model was consulted in ONE batch rather than once per peptide.
        /// </summary>
        private sealed class ScriptedPredictor : IRetentionTimePredictor
        {
            private readonly Dictionary<string, (double? Value, RetentionTimeFailureReason? Reason)> _script;

            internal ScriptedPredictor(Dictionary<string, (double?, RetentionTimeFailureReason?)> script)
                => _script = script;

            internal int SingleCalls { get; private set; }
            internal int BatchCalls { get; private set; }
            internal int PeptidesSeen { get; private set; }

            public string PredictorName => "Scripted";
            public SeparationType SeparationType => SeparationType.HPLC;

            public double? PredictRetentionTimeEquivalent(
                IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
            {
                SingleCalls++;
                if (_script.TryGetValue(peptide.FullSequence, out var scripted))
                {
                    failureReason = scripted.Reason;
                    return scripted.Value;
                }

                failureReason = RetentionTimeFailureReason.PredictionError;
                return null;
            }

            public IReadOnlyList<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)>
                PredictRetentionTimeEquivalents(IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1)
            {
                BatchCalls++;
                var list = peptides.ToList();
                PeptidesSeen += list.Count;
                return list
                    .Select(p => _script.TryGetValue(p.FullSequence, out var s)
                        ? (s.Value, p, s.Reason)
                        : ((double?)null, p, (RetentionTimeFailureReason?)RetentionTimeFailureReason.PredictionError))
                    .ToList();
            }

            public string GetFormattedSequence(
                IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
            {
                failureReason = null;
                return peptide.FullSequence;
            }

            public void Dispose() { }
        }

        private static StubPeptide Peptide(string fullSequence)
            => new() { BaseSequence = fullSequence, FullSequence = fullSequence, FullSequenceWithMassShifts = fullSequence };

        /// <summary>
        /// The contract is FOUR states, not two: a value and a failure reason can arrive TOGETHER. A
        /// predictor in UsePrimarySequence mode reports why it could not represent the peptidoform exactly
        /// and still returns a usable prediction from the primary sequence. Caching only the double would
        /// throw that diagnostic away, which is precisely what PsmData.HasHydrophobicity exists to carry.
        /// </summary>
        [Test]
        public static void WarmedLookupKeepsValueAndFailureReasonTogether()
        {
            var script = new Dictionary<string, (double?, RetentionTimeFailureReason?)>
            {
                ["PLAIN"] = (12.5, null),                                                     // success
                ["UNPREDICTABLE"] = (null, RetentionTimeFailureReason.IncompatibleModifications), // failure
                ["BOTH"] = (34.0, RetentionTimeFailureReason.IncompatibleModifications),      // value AND reason
                ["QUIET"] = (null, null),                                                     // neither
            };
            var inner = new ScriptedPredictor(script);
            var peptides = script.Keys.Select(Peptide).ToList();

            var warm = PrewarmedRetentionTimePredictor.Warm(inner, peptides, maxThreads: 1);

            double? plain = warm.PredictRetentionTimeEquivalent(Peptide("PLAIN"), out var plainReason);
            Assert.That(plain, Is.EqualTo(12.5));
            Assert.That(plainReason, Is.Null);

            double? unpredictable = warm.PredictRetentionTimeEquivalent(Peptide("UNPREDICTABLE"), out var unpredictableReason);
            Assert.That(unpredictable, Is.Null);
            Assert.That(unpredictableReason, Is.EqualTo(RetentionTimeFailureReason.IncompatibleModifications));

            // The one that a value-only cache would corrupt.
            double? both = warm.PredictRetentionTimeEquivalent(Peptide("BOTH"), out var bothReason);
            Assert.That(both, Is.EqualTo(34.0), "a warmed prediction must keep its value");
            Assert.That(bothReason, Is.EqualTo(RetentionTimeFailureReason.IncompatibleModifications),
                "a warmed prediction must keep its failure reason even when it also has a value");

            double? quiet = warm.PredictRetentionTimeEquivalent(Peptide("QUIET"), out var quietReason);
            Assert.That(quiet, Is.Null);
            Assert.That(quietReason, Is.Null);
        }

        /// <summary>
        /// The whole point of the change: the model is consulted once for the batch, not once per peptide,
        /// and repeated peptidoforms are only predicted once.
        /// </summary>
        [Test]
        public static void WarmingPredictsEachDistinctSequenceOnceAndInOneBatch()
        {
            var script = new Dictionary<string, (double?, RetentionTimeFailureReason?)>
            {
                ["AAA"] = (1.0, null),
                ["BBB"] = (2.0, null),
            };
            var inner = new ScriptedPredictor(script);

            // Five peptides, two distinct sequences.
            var peptides = new[] { "AAA", "BBB", "AAA", "AAA", "BBB" }.Select(Peptide).ToList();

            var warm = PrewarmedRetentionTimePredictor.Warm(inner, peptides, maxThreads: 1);

            Assert.That(inner.BatchCalls, Is.EqualTo(1), "one batched call, not one per peptide");
            Assert.That(inner.PeptidesSeen, Is.EqualTo(2), "duplicates must not reach the model");
            Assert.That(inner.SingleCalls, Is.Zero, "the per-peptide path is what this type exists to avoid");
            Assert.That(warm.WarmedSequenceCount, Is.EqualTo(2));

            // And every repeated lookup is served without going back to the model.
            foreach (var p in peptides)
            {
                Assert.That(warm.PredictRetentionTimeEquivalent(p, out _), Is.Not.Null);
            }

            Assert.That(inner.SingleCalls, Is.Zero, "warmed lookups must not fall through");
        }

        /// <summary>
        /// PredictionError is where a blanket catch puts transient native faults, so it must NOT be
        /// memoized -- one bad moment would otherwise poison that peptidoform for the rest of the run.
        /// Every other reason is a pure function of the peptidoform and is safe to keep.
        /// </summary>
        [Test]
        public static void TransientPredictionErrorIsRetriedRatherThanCached()
        {
            var script = new Dictionary<string, (double?, RetentionTimeFailureReason?)>
            {
                ["TRANSIENT"] = (null, RetentionTimeFailureReason.PredictionError),
                ["DETERMINISTIC"] = (null, RetentionTimeFailureReason.SequenceTooLong),
            };
            var inner = new ScriptedPredictor(script);
            var peptides = script.Keys.Select(Peptide).ToList();

            var warm = PrewarmedRetentionTimePredictor.Warm(inner, peptides, maxThreads: 1);

            Assert.That(warm.WarmedSequenceCount, Is.EqualTo(1),
                "only the deterministic failure should have been kept");

            warm.PredictRetentionTimeEquivalent(Peptide("DETERMINISTIC"), out var deterministicReason);
            Assert.That(inner.SingleCalls, Is.Zero, "a deterministic failure is served from the warm set");
            Assert.That(deterministicReason, Is.EqualTo(RetentionTimeFailureReason.SequenceTooLong));

            warm.PredictRetentionTimeEquivalent(Peptide("TRANSIENT"), out var transientReason);
            Assert.That(inner.SingleCalls, Is.EqualTo(1), "a transient failure must be retried against the model");
            Assert.That(transientReason, Is.EqualTo(RetentionTimeFailureReason.PredictionError));
        }

        /// <summary>
        /// Anything not warmed must behave exactly as it did before, so the change cannot alter results for
        /// a peptidoform the warm pass never saw.
        /// </summary>
        [Test]
        public static void UnwarmedPeptideFallsThroughToTheWrappedPredictor()
        {
            var script = new Dictionary<string, (double?, RetentionTimeFailureReason?)>
            {
                ["KNOWN"] = (5.0, null),
                ["LATECOMER"] = (7.0, null),
            };
            var inner = new ScriptedPredictor(script);

            var warm = PrewarmedRetentionTimePredictor.Warm(inner, new[] { Peptide("KNOWN") }, maxThreads: 1);

            Assert.That(warm.PredictRetentionTimeEquivalent(Peptide("LATECOMER"), out var reason), Is.EqualTo(7.0));
            Assert.That(reason, Is.Null);
            Assert.That(inner.SingleCalls, Is.EqualTo(1));
        }

        /// <summary>
        /// FullSequence is annotated non-nullable, but the interface's own documentation says it may be
        /// absent and the shipped implementations defend against that. So must this one.
        /// </summary>
        [Test]
        public static void PeptideWithoutAFullSequenceIsNotWarmedAndDoesNotThrow()
        {
            var inner = new ScriptedPredictor(new Dictionary<string, (double?, RetentionTimeFailureReason?)>());
            var blank = new StubPeptide { BaseSequence = "AAA", FullSequence = string.Empty };

            var warm = PrewarmedRetentionTimePredictor.Warm(inner, new IRetentionPredictable[] { blank }, maxThreads: 1);

            Assert.That(warm.WarmedSequenceCount, Is.Zero);
            Assert.DoesNotThrow(() => warm.PredictRetentionTimeEquivalent(blank, out _));
        }

        /// <summary>
        /// The batched entry point must serve warmed peptidoforms from the table and send ONLY the misses to
        /// the model, in one call, so a later batched consumer cannot reintroduce per-peptide traffic.
        /// </summary>
        [Test]
        public static void BatchedLookupSendsOnlyUnwarmedPeptidesToTheModel()
        {
            var script = new Dictionary<string, (double?, RetentionTimeFailureReason?)>
            {
                ["WARM"] = (3.0, RetentionTimeFailureReason.IncompatibleModifications),
                ["COLD1"] = (4.0, null),
                ["COLD2"] = (5.0, null),
            };
            var inner = new ScriptedPredictor(script);
            var warm = PrewarmedRetentionTimePredictor.Warm(inner, new[] { Peptide("WARM") }, maxThreads: 1);
            int batchCallsAfterWarming = inner.BatchCalls;
            int peptidesSeenAfterWarming = inner.PeptidesSeen;

            var results = warm.PredictRetentionTimeEquivalents(new[] { "WARM", "COLD1", "COLD2" }.Select(Peptide), maxThreads: 1)
                .ToDictionary(r => r.Peptide.FullSequence);

            Assert.That(inner.BatchCalls - batchCallsAfterWarming, Is.EqualTo(1), "the misses go to the model in one call");
            Assert.That(inner.PeptidesSeen - peptidesSeenAfterWarming, Is.EqualTo(2), "a warmed peptidoform must not reach the model");
            Assert.That(inner.SingleCalls, Is.Zero);
            Assert.That(results["WARM"].PredictedValue, Is.EqualTo(3.0));
            Assert.That(results["WARM"].FailureReason, Is.EqualTo(RetentionTimeFailureReason.IncompatibleModifications));
            Assert.That(results["COLD1"].PredictedValue, Is.EqualTo(4.0));
            Assert.That(results["COLD2"].PredictedValue, Is.EqualTo(5.0));

            // Everything warmed: no call reaches the model at all.
            warm.PredictRetentionTimeEquivalents(new[] { Peptide("WARM") });
            Assert.That(inner.BatchCalls - batchCallsAfterWarming, Is.EqualTo(1));
        }

        /// <summary>
        /// The encode phase allocates per peptide before inference starts, so the warm set is handed over in
        /// bounded chunks rather than in one call.
        /// </summary>
        [Test]
        public static void WarmingALargeSetIsChunked()
        {
            var inner = new ScriptedPredictor(new Dictionary<string, (double?, RetentionTimeFailureReason?)>());
            var peptides = Enumerable.Range(0, 50_001).Select(i => Peptide("P" + i));

            PrewarmedRetentionTimePredictor.Warm(inner, peptides, maxThreads: 1);

            Assert.That(inner.BatchCalls, Is.EqualTo(2));
            Assert.That(inner.PeptidesSeen, Is.EqualTo(50_001));
        }

        /// <summary>The decorator must not change what the predictor says about itself.</summary>
        [Test]
        public static void IdentityAndFormattingComeFromTheWrappedPredictor()
        {
            var inner = new ScriptedPredictor(new Dictionary<string, (double?, RetentionTimeFailureReason?)>());
            var warm = PrewarmedRetentionTimePredictor.Warm(inner, Array.Empty<IRetentionPredictable>(), maxThreads: 1);

            Assert.That(warm.PredictorName, Is.EqualTo("Scripted"));
            Assert.That(warm.SeparationType, Is.EqualTo(SeparationType.HPLC));
            Assert.That(warm.GetFormattedSequence(Peptide("AAA"), out var reason), Is.EqualTo("AAA"));
            Assert.That(reason, Is.Null);
        }

        [Test]
        public static void WarmingANullPredictorThrows()
        {
            Assert.Throws<ArgumentNullException>(
                () => PrewarmedRetentionTimePredictor.Warm(null, new[] { Peptide("AAA") }, maxThreads: 1));
        }

        /// <summary>
        /// The wrapped predictor is shared process-wide and is not owned here; disposing it from a borrowed
        /// reference would race its own Dispose, which does not take the model lock.
        /// </summary>
        [Test]
        public static void DisposingTheDecoratorLeavesTheWrappedPredictorUsable()
        {
            var script = new Dictionary<string, (double?, RetentionTimeFailureReason?)> { ["AAA"] = (1.0, null) };
            var inner = new ScriptedPredictor(script);
            var warm = PrewarmedRetentionTimePredictor.Warm(inner, new[] { Peptide("AAA") }, maxThreads: 1);

            warm.Dispose();

            Assert.That(inner.PredictRetentionTimeEquivalent(Peptide("AAA"), out _), Is.EqualTo(1.0));
        }
    }
}
