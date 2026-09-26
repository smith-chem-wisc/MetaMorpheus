using Chromatography;
using Chromatography.RetentionTimePrediction;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using NUnit.Framework;
using Omics;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using TaskLayer;
using UsefulProteomicsDatabases;

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

            /// <summary>Makes the batched call throw, the way a fault inside Chronologer's batched override would.</summary>
            internal bool BatchThrows { get; init; }

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
                if (BatchThrows)
                {
                    throw new AggregateException(new InvalidOperationException("scripted batch failure"));
                }
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

        /// <summary>
        /// Predicts a value for every peptidoform, counts both entry points thread-safely, and records the
        /// most calls it ever saw in flight at once, so a test can prove calls were serialised.
        /// </summary>
        private sealed class CountingPredictor : IRetentionTimePredictor
        {
            private int _singleCalls;
            private int _batchCalls;
            private int _peptidesSeen;
            private int _inFlight;
            private int _maxConcurrentCalls;

            internal int SingleCallDelayMs { get; init; }
            internal int SingleCalls => _singleCalls;
            internal int BatchCalls => _batchCalls;
            internal int PeptidesSeen => _peptidesSeen;
            internal int MaxConcurrentCalls => _maxConcurrentCalls;

            public string PredictorName => "Counting";
            public SeparationType SeparationType => SeparationType.HPLC;

            public double? PredictRetentionTimeEquivalent(
                IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
            {
                Interlocked.Increment(ref _singleCalls);
                int inFlight = Interlocked.Increment(ref _inFlight);
                int seen;
                while (inFlight > (seen = Volatile.Read(ref _maxConcurrentCalls))
                       && Interlocked.CompareExchange(ref _maxConcurrentCalls, inFlight, seen) != seen)
                {
                }
                if (SingleCallDelayMs > 0)
                {
                    Thread.Sleep(SingleCallDelayMs);
                }
                Interlocked.Decrement(ref _inFlight);

                failureReason = null;
                return Predict(peptide);
            }

            public IReadOnlyList<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)>
                PredictRetentionTimeEquivalents(IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1)
            {
                Interlocked.Increment(ref _batchCalls);
                var list = peptides.ToList();
                Interlocked.Add(ref _peptidesSeen, list.Count);
                return list.Select(p => (Predict(p), p, (RetentionTimeFailureReason?)null)).ToList();
            }

            private static double? Predict(IRetentionPredictable peptide) => peptide.FullSequence.Length;

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
        /// A warmed lookup returns the pair the BATCHED call gave, in all four combinations. That is not a
        /// promise that it matches the single-peptide path: Chronologer's batched override returns a null
        /// reason for every success, so a UsePrimarySequence reason the single path would keep never reaches
        /// this table. Nothing reads the reason today; both consumers only ask whether a value exists.
        /// </summary>
        [Test]
        public static void WarmedLookupReturnsThePairTheBatchedCallGave()
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

            double? both = warm.PredictRetentionTimeEquivalent(Peptide("BOTH"), out var bothReason);
            Assert.That(both, Is.EqualTo(34.0));
            Assert.That(bothReason, Is.EqualTo(RetentionTimeFailureReason.IncompatibleModifications),
                "the table stores what the batched call returned, value and reason together");

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
        /// Koina marks every row of a failed HTTP batch PredictionError, so a batched PredictionError is not
        /// kept from warming: it gets ONE retry on the single path, and that answer is then remembered, so
        /// the peptidoform does not go back to the model once per hypothesis.
        /// </summary>
        [Test]
        public static void BatchedPredictionErrorIsRetriedOnceThenRemembered()
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
            Assert.That(inner.SingleCalls, Is.EqualTo(1), "a batched PredictionError must be retried against the model");
            Assert.That(transientReason, Is.EqualTo(RetentionTimeFailureReason.PredictionError));

            warm.PredictRetentionTimeEquivalent(Peptide("TRANSIENT"), out var secondReason);
            Assert.That(inner.SingleCalls, Is.EqualTo(1), "the retry's answer is remembered, not asked again");
            Assert.That(secondReason, Is.EqualTo(RetentionTimeFailureReason.PredictionError));
            Assert.That(warm.MissCount, Is.EqualTo(1));
        }

        /// <summary>
        /// Chronologer's batched override has no catch, unlike the single path, which turns an exception into
        /// one PredictionError. A batch that throws must cost that batch's warming, not the run.
        /// </summary>
        [Test]
        public static void ABatchThatThrowsIsLeftUnwarmedAndFallsThrough()
        {
            var script = new Dictionary<string, (double?, RetentionTimeFailureReason?)> { ["AAA"] = (1.0, null) };
            var inner = new ScriptedPredictor(script) { BatchThrows = true };

            PrewarmedRetentionTimePredictor warm = null;
            Assert.DoesNotThrow(() => warm = PrewarmedRetentionTimePredictor.Warm(inner, new[] { Peptide("AAA") }, maxThreads: 1));

            Assert.That(warm.WarmedSequenceCount, Is.Zero);
            Assert.That(warm.FailedChunkCount, Is.EqualTo(1));
            Assert.That(warm.FirstChunkFailure, Is.EqualTo("scripted batch failure"));
            Assert.That(warm.PredictRetentionTimeEquivalent(Peptide("AAA"), out _), Is.EqualTo(1.0));
            Assert.That(warm.PredictRetentionTimeEquivalent(Peptide("AAA"), out _), Is.EqualTo(1.0));
            Assert.That(inner.SingleCalls, Is.EqualTo(1), "the fall-through is asked once per peptidoform");
        }

        /// <summary>
        /// Misses arrive from PEP's Parallel.ForEach, and Koina's RetentionTimeModel is documented as not
        /// thread-safe. Each missed peptidoform must reach the wrapped predictor once, and never concurrently.
        /// </summary>
        [Test]
        public static void ConcurrentMissesReachTheWrappedPredictorOnceEachAndNeverAtOnce()
        {
            var inner = new CountingPredictor { SingleCallDelayMs = 1 };
            var warm = PrewarmedRetentionTimePredictor.Warm(inner, Array.Empty<IRetentionPredictable>(), maxThreads: 1);
            var sequences = Enumerable.Range(0, 20).Select(i => "SEQ" + i).ToArray();

            Parallel.For(0, 2_000, new ParallelOptions { MaxDegreeOfParallelism = 16 },
                i => warm.PredictRetentionTimeEquivalent(Peptide(sequences[i % sequences.Length]), out _));

            Assert.That(inner.SingleCalls, Is.EqualTo(sequences.Length));
            Assert.That(warm.MissCount, Is.EqualTo(sequences.Length));
            Assert.That(inner.MaxConcurrentCalls, Is.EqualTo(1));
        }

        /// <summary>
        /// The regression the log line exists to reveal, pinned end to end: after warming, a full PEP pass
        /// never asks the model about one peptide at a time.
        /// </summary>
        [Test]
        public static void PepAnalysisMakesNoSinglePeptideCallsAfterWarming()
        {
            var (psms, fsp) = SearchForPep("HPLC");
            var inner = new CountingPredictor();

            string log = new PepAnalysisEngine(psms, "standard", fsp, PepOutputFolder(), inner).ComputePEPValuesForAllPSMs();

            Assert.That(inner.BatchCalls, Is.EqualTo(1));
            Assert.That(inner.SingleCalls, Is.Zero, "every retention time must come from the warm set");
            Assert.That(log, Does.Contain(" distinct peptidoforms, in batches, and 0 one at a time"));
        }

        /// <summary>
        /// A CZE file is scored on mobility, so the only retention times read for it are the q &lt;= 0.01
        /// targets in the reference distribution. Warm exactly those, and still make no single-peptide calls.
        /// </summary>
        [Test]
        public static void CzeFilesWarmOnlyTheReferenceTargets()
        {
            var (psms, fsp) = SearchForPep("CZE");
            var inner = new CountingPredictor();
            int referenceTargets = psms
                .Where(p => !p.IsDecoy && p.FdrInfo.QValue <= 0.01)
                .SelectMany(p => p.BestMatchingBioPolymersWithSetMods)
                .Select(h => h.SpecificBioPolymer)
                .OfType<PeptideWithSetModifications>()
                .Select(p => p.FullSequence)
                .Distinct()
                .Count();
            int everyHypothesis = psms
                .SelectMany(p => p.BestMatchingBioPolymersWithSetMods)
                .Select(h => h.SpecificBioPolymer.FullSequence)
                .Distinct()
                .Count();
            Assert.That(referenceTargets, Is.GreaterThan(0).And.LessThan(everyHypothesis), "the data must contain PSMs the filter drops");

            new PepAnalysisEngine(psms, "standard", fsp, PepOutputFolder(), inner).ComputePEPValuesForAllPSMs();

            Assert.That(inner.PeptidesSeen, Is.EqualTo(referenceTargets));
            Assert.That(inner.SingleCalls, Is.Zero);
        }

        /// <summary>
        /// ScriptedPredictor cannot show that the shipped predictor's batched answers match its single-peptide
        /// answers. Warm the real Chronologer and compare, on modified and unmodified peptidoforms.
        /// </summary>
        [Test]
        public static void WarmedChronologerAgreesWithItsSinglePeptidePath()
        {
            IRetentionTimePredictor chronologer = FdrAnalysisEngine.GetChronologer();
            Assert.That(chronologer, Is.Not.Null, "Chronologer's native libraries could not be loaded");

            Modification oxidation = GlobalVariables.AllModsKnown.First(m => m.IdWithMotif == "Oxidation on M");
            Modification carbamidomethyl = GlobalVariables.AllModsKnown.First(m => m.IdWithMotif == "Carbamidomethyl on C");
            var protein = new Protein("MPEPTIDEKCAMSTERDAMKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKMLLVGGAR", "P0");
            List<IRetentionPredictable> peptides = protein
                .Digest(new DigestionParams(), new List<Modification> { carbamidomethyl }, new List<Modification> { oxidation })
                .Cast<IRetentionPredictable>()
                .ToList();
            Assert.That(peptides.Any(p => p.FullSequence.Contains("Oxidation")), "need a variably modified peptidoform");

            var warm = PrewarmedRetentionTimePredictor.Warm(chronologer, peptides, maxThreads: 2);

            int compared = 0;
            foreach (IRetentionPredictable peptide in peptides)
            {
                double? single = chronologer.PredictRetentionTimeEquivalent(peptide, out _);
                double? warmed = warm.PredictRetentionTimeEquivalent(peptide, out _);
                if (single.HasValue)
                {
                    Assert.That(warmed, Is.EqualTo(single.Value).Within(1e-4), peptide.FullSequence);
                    compared++;
                }
                else
                {
                    Assert.That(warmed, Is.Null, peptide.FullSequence);
                }
            }

            Assert.That(compared, Is.GreaterThan(0));
            Assert.That(warm.FailedChunkCount, Is.Zero);
            Assert.That(warm.MissCount, Is.Zero);
        }

        private static string PepOutputFolder() => Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\");

        /// <summary>The same search <see cref="FdrTest.TestComputePEPValue"/> runs, with the file's separation type set.</summary>
        private static (List<SpectralMatch> Psms, List<(string fileName, CommonParameters fileSpecificParameters)> Fsp) SearchForPep(string separationType)
        {
            string origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_HeLa_04_subset_longestSeq.mzML");
            var commonParameters = new CommonParameters(digestionParams: new DigestionParams(), separationType: separationType);
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>
            {
                ("TaGe_SA_HeLa_04_subset_longestSeq.mzML", commonParameters)
            };

            var myMsDataFile = new MyFileManager(true).LoadFile(origDataFile, commonParameters);
            List<Protein> proteinList = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\hela_snip_for_unitTest.fasta"), true, DecoyType.Reverse, false, out _,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex, -1);
            var scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, origDataFile, commonParameters).OrderBy(b => b.PrecursorMass).ToArray();
            SpectralMatch[] allPsms = new PeptideSpectralMatch[scans.Length];
            new ClassicSearchEngine(allPsms, scans, new List<Modification>(), new List<Modification>(), null, null, null,
                proteinList, new SinglePpmAroundZeroSearchMode(5), commonParameters, fsp, null, new List<string>(), false).Run();
            List<SpectralMatch> psms = allPsms.Where(p => p != null).ToList();
            new FdrAnalysisEngine(psms, 1, commonParameters, fsp, new List<string>(), doPEP: false).Run();
            return (psms, fsp);
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
