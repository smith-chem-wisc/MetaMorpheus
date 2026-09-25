using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using Microsoft.ML;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// Iterative (semi-supervised) PEP training. Every test searches the same small HeLa subset afresh, because
    /// the PEP engine writes onto the matches it is given.
    /// </summary>
    [TestFixture]
    public static class PepIterativeTrainingTests
    {
        private const string DataFileName = "TaGe_SA_HeLa_04_subset_longestSeq.mzML";

        private static string OutputFolder => Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPepIterativeTraining");

        private static List<(string fileName, CommonParameters fileSpecificParameters)> FileSpecificParameters(CommonParameters commonParameters)
            => new() { (DataFileName, commonParameters) };

        private static List<SpectralMatch> SearchHelaSubset(out CommonParameters commonParameters)
        {
            commonParameters = new CommonParameters(digestionParams: new DigestionParams());
            var fsp = FileSpecificParameters(commonParameters);
            var dataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData", DataFileName);
            var msDataFile = new MyFileManager(true).LoadFile(dataFile, commonParameters);
            List<Protein> proteins = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\hela_snip_for_unitTest.fasta"),
                true, DecoyType.Reverse, false, out _, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex,
                ProteinDbLoader.UniprotGeneNameRegex, ProteinDbLoader.UniprotOrganismRegex, -1);
            var scans = MetaMorpheusTask.GetMs2Scans(msDataFile, dataFile, commonParameters).OrderBy(b => b.PrecursorMass).ToArray();
            SpectralMatch[] psmArray = new PeptideSpectralMatch[scans.Length];
            new ClassicSearchEngine(psmArray, scans, new List<Modification>(), new List<Modification>(), null, null, null,
                proteins, new SinglePpmAroundZeroSearchMode(5), commonParameters, fsp, null, new List<string>(), false).Run();
            var psms = psmArray.Where(p => p != null).ToList();
            // q-values only: fewer than 1000 matches, so FdrAnalysisEngine does not run PEP itself
            new FdrAnalysisEngine(psms, 1, commonParameters, fsp, new List<string>()).Run();
            return psms;
        }

        private static PepAnalysisEngine NewEngine(List<SpectralMatch> psms, CommonParameters commonParameters, int maxTrainingRounds)
        {
            Directory.CreateDirectory(OutputFolder);
            return new PepAnalysisEngine(psms, "standard", FileSpecificParameters(commonParameters), OutputFolder)
            {
                MaxTrainingRounds = maxTrainingRounds
            };
        }

        private static double[] Peps(List<SpectralMatch> psms) => psms.Select(p => p.PsmFdrInfo.PEP).ToArray();

        private static int RoundsRun(string metrics)
            => int.Parse(Regex.Match(metrics, @"Training Rounds Run:\s+(\d+)").Groups[1].Value);

        [TearDown]
        public static void TearDown()
        {
            if (Directory.Exists(OutputFolder))
            {
                Directory.Delete(OutputFolder, true);
            }
        }

        /// <summary>
        /// With iteration off, the engine must be the algorithm it was before iteration existed. The reference
        /// below IS that algorithm -- the train-once, four-fold loop -- spelled out from the engine's public
        /// pieces. Compared bit for bit on the same machine, so it does not depend on platform floating point.
        /// </summary>
        [Test]
        public static void OneRound_IsBitIdenticalToTheTrainOnceAlgorithm()
        {
            var psms = SearchHelaSubset(out var commonParameters);
            string metrics = NewEngine(psms, commonParameters, 1).ComputePEPValuesForAllPSMs();
            Assert.That(RoundsRun(metrics), Is.EqualTo(1));

            var referencePsms = SearchHelaSubset(out commonParameters);
            var reference = NewEngine(referencePsms, commonParameters, 1);
            var groups = reference.UsePeptideLevelQValueForTraining
                ? SpectralMatchGroup.GroupByBaseSequence(reference.AllPsms)
                : SpectralMatchGroup.GroupByIndividualPsm(reference.AllPsms);
            if (reference.UsePeptideLevelQValueForTraining && (groups.Count(g => g.BestMatch.IsDecoy) < 4 || groups.Count(g => !g.BestMatch.IsDecoy) < 4))
            {
                groups = SpectralMatchGroup.GroupByIndividualPsm(reference.AllPsms);
                reference.UsePeptideLevelQValueForTraining = false;
            }
            var indices = PepAnalysisEngine.GetPeptideGroupIndices(groups, 4);
            var data = Enumerable.Range(0, 4).Select(g => reference.CreatePsmData("standard", groups, indices[g])).ToArray();
            var mlContext = new MLContext(seed: 42);
            var pipeline = mlContext.Transforms.Concatenate("Features", reference.TrainingVariables)
                .Append(mlContext.BinaryClassification.Trainers.FastTree(reference.BGDTreeOptions));
            for (int fold = 0; fold < 4; fold++)
            {
                var others = Enumerable.Range(0, 4).Where(g => g != fold).ToList();
                var model = pipeline.Fit(mlContext.Data.LoadFromEnumerable(data[others[0]].Concat(data[others[1]].Concat(data[others[2]]))));
                mlContext.BinaryClassification.Evaluate(model.Transform(mlContext.Data.LoadFromEnumerable(data[fold])), "Label", "Score");
                reference.Compute_PSM_PEP(groups, indices[fold], mlContext, model, "standard", OutputFolder);
            }

            Assert.That(Peps(psms), Is.EqualTo(Peps(referencePsms)));
        }

        /// <summary>
        /// Whatever round the loop stops at, and whether it stops by the tolerance or by rejecting a worse
        /// round, the PEPs it leaves must be exactly the last KEPT round's. So a run capped at the number of
        /// rounds an uncapped run kept must give identical PEPs. When the uncapped run ended by rejecting a
        /// round, this is the restore-on-no-improvement guard at work.
        /// </summary>
        [Test]
        public static void Iterating_LeavesExactlyTheLastKeptRoundsPeps()
        {
            var psms = SearchHelaSubset(out var commonParameters);
            string metrics = NewEngine(psms, commonParameters, PepAnalysisEngine.IterativeTrainingRoundCap).ComputePEPValuesForAllPSMs();
            int kept = RoundsRun(metrics);
            string progress = File.ReadAllText(Path.Combine(OutputFolder, "pep_training_rounds.txt"));
            TestContext.WriteLine(progress);
            // The last line for a round is either the rejected round after the kept ones, or the kept round that stopped.
            bool reverted = progress.Contains($"round {kept}: accepted") && progress.Contains(nameof(PepAnalysisEngine.RoundVerdict.RevertAndStop));
            bool stopped = progress.Contains($"round {kept - 1}: accepted") && progress.Contains(nameof(PepAnalysisEngine.RoundVerdict.KeepAndStop));
            Assert.That(reverted || stopped, progress);

            var cappedPsms = SearchHelaSubset(out commonParameters);
            string cappedMetrics = NewEngine(cappedPsms, commonParameters, kept).ComputePEPValuesForAllPSMs();

            Assert.That(RoundsRun(cappedMetrics), Is.EqualTo(kept));
            Assert.That(Peps(psms), Is.EqualTo(Peps(cappedPsms)));
            Assert.That(Peps(psms).All(p => p >= 0 && p <= 1));
        }

        /// <summary>
        /// A later round in which a fold has no positive examples must keep the previous round's PEPs
        /// untouched: the check runs for every fold before any model is refitted or any PEP is overwritten.
        /// </summary>
        [Test]
        public static void LaterRoundWithoutPositives_KeepsThePreviousRound()
        {
            var oneRoundPsms = SearchHelaSubset(out var commonParameters);
            NewEngine(oneRoundPsms, commonParameters, 1).ComputePEPValuesForAllPSMs();

            var psms = SearchHelaSubset(out commonParameters);
            var engine = NewEngine(psms, commonParameters, PepAnalysisEngine.IterativeTrainingRoundCap);
            engine.SelectNextRoundPositives = (matches, pep) => new HashSet<SpectralMatch>();
            string metrics = engine.ComputePEPValuesForAllPSMs();

            Assert.That(RoundsRun(metrics), Is.EqualTo(1));
            Assert.That(Peps(psms), Is.EqualTo(Peps(oneRoundPsms)));
            Assert.That(File.ReadAllText(Path.Combine(OutputFolder, "pep_training_rounds.txt")), Does.Contain("keeping round 0"));
        }

        /// <summary>
        /// One group with no positive examples while the other three have some. Every fold's training set then
        /// has both classes, but one held-out fold does not. The engine must refuse cleanly before anything
        /// trains, and leave no PEP assigned.
        /// </summary>
        [Test]
        public static void AGroupWithoutPositives_FailsCleanly_BeforeAnyPepIsAssigned()
        {
            var psms = SearchHelaSubset(out var commonParameters);
            var engine = NewEngine(psms, commonParameters, PepAnalysisEngine.IterativeTrainingRoundCap);

            // Keep positives in groups 1-3 only. Group membership is the engine's own partition.
            var groups = engine.UsePeptideLevelQValueForTraining
                ? SpectralMatchGroup.GroupByBaseSequence(engine.AllPsms)
                : SpectralMatchGroup.GroupByIndividualPsm(engine.AllPsms);
            var indices = PepAnalysisEngine.GetPeptideGroupIndices(groups, 4);
            foreach (var psm in indices[0].SelectMany(i => groups[i]).Where(p => !p.IsDecoy))
            {
                psm.PsmFdrInfo.QValue = 1;
                psm.PeptideFdrInfo.QValue = 1;
            }
            Assume.That(indices.Skip(1).All(g => g.SelectMany(i => groups[i].GetBestMatches())
                .Any(p => !p.IsDecoy && p.GetFdrInfo(engine.UsePeptideLevelQValueForTraining).QValue <= engine.QValueCutoff)));
            double[] before = Peps(psms);

            string result = engine.ComputePEPValuesForAllPSMs();

            Assert.That(result, Is.EqualTo("Posterior error probability analysis failed. This can occur for small data sets when some sample groups are missing positive or negative training examples."));
            Assert.That(Peps(psms), Is.EqualTo(before));
        }

        /// <summary>
        /// The stopping count must be the number FdrAnalysisEngine reports: target peptides whose PEP q-value
        /// is below the cutoff, with the PEP-best match per full sequence.
        /// </summary>
        [Test]
        public static void CountAcceptedPeptides_MatchesTheReportedPeptidePepQValues()
        {
            var property = typeof(FdrAnalysisEngine).GetProperty("QvalueThresholdOverride")!;
            try
            {
                property.SetValue(null, true);
                var psms = SearchHelaSubset(out var commonParameters);
                new FdrAnalysisEngine(psms, 1, commonParameters, FileSpecificParameters(commonParameters), new List<string>()).Run();

                var peptides = psms
                    .OrderBy(p => p.FdrInfo.PEP)
                    .ThenByDescending(p => p)
                    .GroupBy(p => p.FullSequence)
                    .Select(g => g.First())
                    .ToList();

                // This subset is too small for (decoys + 1) / targets to reach 0.01, so compare across cutoffs.
                int acceptedSomewhere = 0;
                foreach (double cutoff in new[] { 0.01, 0.05, 0.1, 0.2, 0.5 })
                {
                    int reported = peptides.Count(p => !p.IsDecoy && p.PeptideFdrInfo.PEP_QValue < cutoff);
                    Assert.That(PepAnalysisEngine.CountAcceptedPeptides(psms, cutoff), Is.EqualTo(reported), $"cutoff {cutoff}");
                    acceptedSomewhere = Math.Max(acceptedSomewhere, reported);
                }
                Assert.That(acceptedSomewhere, Is.GreaterThan(0));
            }
            finally
            {
                property.SetValue(null, false);
            }
        }

        /// <summary>
        /// The next round's positives are cut by RANK. Every match below is tied at the same PEP, so a
        /// threshold rule (PEP &lt;= t) would admit every target. The rank cut admits only the targets the sweep
        /// ordered ahead of the first decoy, which is where q first exceeds a cutoff of 0.
        /// </summary>
        [Test]
        public static void SelectPositives_CutsTiesByRank()
        {
            var ordered = SearchHelaSubset(out _).OrderByDescending(p => p).ToList();
            int decoyIndex = Enumerable.Range(0, ordered.Count).First(i => ordered[i].IsDecoy
                && ordered.Take(i).Count(p => !p.IsDecoy) >= 3 && ordered.Skip(i + 1).Any(p => !p.IsDecoy));
            var matches = ordered.Take(decoyIndex).Where(p => !p.IsDecoy).TakeLast(3)
                .Append(ordered[decoyIndex])
                .Concat(ordered.Skip(decoyIndex + 1).Where(p => !p.IsDecoy).Take(3))
                .Reverse() // input order must not matter
                .ToList();
            var tiedPep = new double[matches.Count];

            var positives = PepAnalysisEngine.SelectPositives(matches, tiedPep, 0.0);

            var expected = matches.OrderByDescending(m => m).TakeWhile(m => !m.IsDecoy).ToList();
            Assert.That(expected.Count, Is.EqualTo(3));
            Assert.That(matches.Count(m => !m.IsDecoy), Is.GreaterThan(3), "a PEP <= t rule would have admitted these too");
            Assert.That(positives, Is.EquivalentTo(expected));
        }

        [Test]
        [TestCase(0, 100, 0, 1, "KeepAndStop", TestName = "OneRoundCap_StopsAfterRoundZero")]
        [TestCase(0, 100, 0, 10, "Continue", TestName = "RoundZero_Continues")]
        [TestCase(1, 99, 100, 10, "RevertAndStop", TestName = "WorseRound_IsReverted")]
        [TestCase(1, 100, 100, 10, "RevertAndStop", TestName = "EqualRound_IsReverted")]
        [TestCase(1, 10000, 9999, 10, "KeepAndStop", TestName = "GainBelowTolerance_KeepsAndStops")]
        [TestCase(1, 1010, 1000, 10, "Continue", TestName = "GainAtLeastTolerance_Continues")]
        [TestCase(9, 2000, 1000, 10, "KeepAndStop", TestName = "Cap_KeepsAndStops")]
        [TestCase(1, 5, 0, 10, "Continue", TestName = "GainFromZero_Continues")]
        public static void JudgeRound(int round, int accepted, int previousAccepted, int maxRounds, string expected)
        {
            Assert.That(PepAnalysisEngine.JudgeRound(round, accepted, previousAccepted, 0.001, maxRounds).ToString(), Is.EqualTo(expected));
        }

        [Test]
        public static void IterativePepTraining_IsOffByDefault()
        {
            Assert.That(new SearchParameters().IterativePepTraining, Is.False);
            Assert.That(new PepAnalysisEngine(SearchHelaSubset(out var commonParameters), "standard", FileSpecificParameters(commonParameters), null).MaxTrainingRounds, Is.EqualTo(1));
        }
    }
}
