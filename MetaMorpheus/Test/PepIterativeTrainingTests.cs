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
    [NonParallelizable] // sets FdrAnalysisEngine.QvalueThresholdOverride, a process-wide static
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

        /// <summary>
        /// The round count in the results block. A train-once run carries no round lines at all, so that is 1.
        /// </summary>
        private static int RoundsRun(string metrics)
        {
            Match match = Regex.Match(metrics, @"Training Rounds Run:\s+(\d+)");
            return match.Success ? int.Parse(match.Groups[1].Value) : 1;
        }

        private static string ProgressFile => Path.Combine(OutputFolder, "pep_training_rounds.txt");

        /// <summary>
        /// Scripts the stopping count, one value per round, so a test can decide which rounds are kept.
        /// </summary>
        private static Func<List<SpectralMatch>, int> ScriptedCounts(params int[] counts)
        {
            int round = 0;
            return _ => counts[round++];
        }

        /// <summary>
        /// Every target as a positive: a label set unlike round 0's search-score cut, so a round trained on it
        /// assigns different PEPs.
        /// </summary>
        private static HashSet<SpectralMatch> AllTargets(List<SpectralMatch> matches, double[] pep)
            => matches.Where(m => !m.IsDecoy).ToHashSet();

        private static PepAnalysisEngine ScriptedEngine(List<SpectralMatch> psms, CommonParameters commonParameters, int maxTrainingRounds,
            bool prune, params int[] counts)
        {
            Directory.CreateDirectory(OutputFolder);
            return new PepAnalysisEngine(psms, "standard", FileSpecificParameters(commonParameters), OutputFolder, pruneAmbiguousHypotheses: prune)
            {
                MaxTrainingRounds = maxTrainingRounds,
                SelectNextRoundPositives = AllTargets,
                CountAcceptedTargets = ScriptedCounts(counts)
            };
        }

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
            // Its output is the train-once engine's too: no round lines in the results block, and no progress file.
            Assert.That(metrics, Does.Not.Contain("Training Rounds Run"));
            Assert.That(metrics, Does.Not.Contain("Accepted After Training"));
            Assert.That(File.Exists(ProgressFile), Is.False);

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
            // The last line for a round is the rejected round after the kept ones, the kept round that stopped, or a round
            // whose relabelling starved a fold (on this small file no target reaches the training cutoff after round 0).
            bool reverted = progress.Contains($"round {kept}: accepted") && progress.Contains(nameof(PepAnalysisEngine.RoundVerdict.RevertAndStop));
            bool stopped = progress.Contains($"round {kept - 1}: accepted") && progress.Contains(nameof(PepAnalysisEngine.RoundVerdict.KeepAndStop));
            bool starved = progress.Contains($"keeping round {kept - 1}");
            Assert.That(reverted || stopped || starved, progress);

            var cappedPsms = SearchHelaSubset(out commonParameters);
            string cappedMetrics = NewEngine(cappedPsms, commonParameters, kept).ComputePEPValuesForAllPSMs();

            Assert.That(RoundsRun(cappedMetrics), Is.EqualTo(kept));
            Assert.That(Peps(psms), Is.EqualTo(Peps(cappedPsms)));
            Assert.That(Peps(psms).All(p => p >= 0 && p <= 1));
        }

        /// <summary>
        /// The paths the small fixture cannot reach on its own, where a round after round 0 is KEPT. With the
        /// stopping count scripted to 10, 20, 15: round 1 is kept, round 2 is rejected, and the PEPs left must be
        /// round 1's. Those are pinned by a run capped at 2 rounds, which keeps round 1 and stops at the cap, and
        /// they must differ from round 0's, or the comparison would say nothing.
        /// </summary>
        [Test]
        public static void LaterKeptRound_IsWhatARejectedRoundRestores()
        {
            var roundZeroPsms = SearchHelaSubset(out var commonParameters);
            NewEngine(roundZeroPsms, commonParameters, 1).ComputePEPValuesForAllPSMs();

            var psms = SearchHelaSubset(out commonParameters);
            string metrics = ScriptedEngine(psms, commonParameters, PepAnalysisEngine.IterativeTrainingRoundCap, false, 10, 20, 15)
                .ComputePEPValuesForAllPSMs();
            string progress = File.ReadAllText(ProgressFile);
            Assert.That(RoundsRun(metrics), Is.EqualTo(2), progress);
            Assert.That(progress, Does.Contain($"round 2: accepted 15  {PepAnalysisEngine.RoundVerdict.RevertAndStop}"));

            var cappedPsms = SearchHelaSubset(out commonParameters);
            string cappedMetrics = ScriptedEngine(cappedPsms, commonParameters, 2, false, 10, 20).ComputePEPValuesForAllPSMs();
            Assert.That(RoundsRun(cappedMetrics), Is.EqualTo(2));

            Assert.That(Peps(psms), Is.EqualTo(Peps(cappedPsms)));
            Assert.That(Peps(psms), Is.Not.EqualTo(Peps(roundZeroPsms)), "round 1 must have changed the PEPs");
        }

        /// <summary>
        /// Iteration and pruning are independent. Pruning waits until the loop has chosen a round, so the PEPs of a
        /// pruning run equal a non-pruning run's round for round, and what it prunes is judged on the kept round.
        /// </summary>
        [Test]
        public static void Pruning_DoesNotLimitIteration_AndPrunesFromTheKeptRound()
        {
            var plainPsms = SearchHelaSubset(out var commonParameters);
            ScriptedEngine(plainPsms, commonParameters, PepAnalysisEngine.IterativeTrainingRoundCap, false, 10, 20, 15)
                .ComputePEPValuesForAllPSMs();
            int hypothesesBefore = plainPsms.Sum(p => p.BestMatchingBioPolymersWithSetMods.Count());

            var prunedPsms = SearchHelaSubset(out commonParameters);
            string prunedMetrics = ScriptedEngine(prunedPsms, commonParameters, PepAnalysisEngine.IterativeTrainingRoundCap, true, 10, 20, 15)
                .ComputePEPValuesForAllPSMs();
            int removed = int.Parse(Regex.Match(prunedMetrics, @"Removed:\s+(\d+)").Groups[1].Value);

            Assert.That(RoundsRun(prunedMetrics), Is.EqualTo(2));
            Assert.That(Peps(prunedPsms), Is.EqualTo(Peps(plainPsms)));
            Assert.That(prunedPsms.Sum(p => p.BestMatchingBioPolymersWithSetMods.Count()), Is.EqualTo(hypothesesBefore - removed));
        }

        /// <summary>
        /// The setting reaches the engine: FdrAnalysisEngine.Compute_PEPValue iterates only when asked to, which
        /// shows as the progress file (written only when iterating).
        /// </summary>
        [Test]
        [TestCase(false)]
        [TestCase(true)]
        public static void ComputePepValue_IteratesOnlyWhenAsked(bool iterativePepTraining)
        {
            var psms = SearchHelaSubset(out var commonParameters);
            var fsp = FileSpecificParameters(commonParameters);
            Directory.CreateDirectory(OutputFolder);
            var results = new FdrAnalysisResults(new FdrAnalysisEngine(psms, 1, commonParameters, fsp, new List<string>()), "PSM");

            FdrAnalysisEngine.Compute_PEPValue(results, psms, fsp, OutputFolder, iterativePepTraining: iterativePepTraining);

            Assert.That(File.Exists(ProgressFile), Is.EqualTo(iterativePepTraining));
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
        /// One positive example per group. That passes the all-groups check, and round 0 trains, but a later round's
        /// labels can leave a held-out fold with no positives. Held-out folds are evaluated against the round-0 labels,
        /// so iterating must finish normally with a PEP on every match instead of throwing in ML.NET's Evaluate.
        /// </summary>
        [Test]
        public static void OnePositivePerGroup_IteratesWithoutThrowing()
        {
            var psms = SearchHelaSubset(out var commonParameters);
            var engine = NewEngine(psms, commonParameters, PepAnalysisEngine.IterativeTrainingRoundCap);

            var groups = engine.UsePeptideLevelQValueForTraining
                ? SpectralMatchGroup.GroupByBaseSequence(engine.AllPsms)
                : SpectralMatchGroup.GroupByIndividualPsm(engine.AllPsms);
            var indices = PepAnalysisEngine.GetPeptideGroupIndices(groups, 4);
            foreach (var group in indices)
            {
                int keep = group.First(i => groups[i].GetBestMatches()
                    .Any(p => !p.IsDecoy && p.GetFdrInfo(engine.UsePeptideLevelQValueForTraining).QValue <= engine.QValueCutoff));
                foreach (var psm in group.Where(i => i != keep).SelectMany(i => groups[i]).Where(p => !p.IsDecoy))
                {
                    psm.PsmFdrInfo.QValue = 1;
                    psm.PeptideFdrInfo.QValue = 1;
                }
            }

            string metrics = engine.ComputePEPValuesForAllPSMs();

            Assert.That(RoundsRun(metrics), Is.GreaterThanOrEqualTo(1));
            Assert.That(Peps(psms).All(p => p >= 0 && p <= 1));
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
        /// threshold rule (PEP &lt;= t) would admit every target. With N targets, one decoy, then M &lt; N targets,
        /// q is (decoys + 1) / targets: 1/N at the last target before the decoy, and at least 2/(N+M) at every
        /// target after it. A cutoff between the two admits exactly the targets the sweep ordered ahead of the decoy.
        /// </summary>
        [Test]
        public static void SelectPositives_CutsTiesByRank()
        {
            // No match mixing target and decoy hypotheses, so each counts as exactly one target or one decoy.
            var ordered = SearchHelaSubset(out _).Where(p => p.BestMatchingBioPolymersWithSetMods.All(h => h.IsDecoy) || p.BestMatchingBioPolymersWithSetMods.All(h => !h.IsDecoy))
                .OrderByDescending(p => p).ToList();
            int decoyIndex = Enumerable.Range(0, ordered.Count).First(i => ordered[i].IsDecoy
                && ordered.Take(i).Count(p => !p.IsDecoy) >= 2 && ordered.Skip(i + 1).Any(p => !p.IsDecoy));
            var before = ordered.Take(decoyIndex).Where(p => !p.IsDecoy).TakeLast(10).ToList();
            var after = ordered.Skip(decoyIndex + 1).Where(p => !p.IsDecoy).Take(Math.Min(3, before.Count - 1)).ToList();
            var matches = before
                .Append(ordered[decoyIndex])
                .Concat(after)
                .Reverse() // input order must not matter
                .ToList();
            var tiedPep = new double[matches.Count];
            double cutoff = (1.0 / before.Count + 2.0 / (before.Count + after.Count)) / 2;

            var positives = PepAnalysisEngine.SelectPositives(matches, tiedPep, cutoff);

            var expected = matches.OrderByDescending(m => m).TakeWhile(m => !m.IsDecoy).ToList();
            Assert.That(expected, Is.EquivalentTo(before));
            Assert.That(matches.Count(m => !m.IsDecoy), Is.GreaterThan(before.Count), "a PEP <= t rule would have admitted these too");
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
        public static void IterativePepTraining_IsOnByDefault_ForASearchTask()
        {
            Assert.That(new SearchParameters().IterativePepTraining, Is.True);
            // The engine on its own trains once; only a search task's setting turns iteration on.
            Assert.That(new PepAnalysisEngine(SearchHelaSubset(out var commonParameters), "standard", FileSpecificParameters(commonParameters), null).MaxTrainingRounds, Is.EqualTo(1));
        }
    }
}
