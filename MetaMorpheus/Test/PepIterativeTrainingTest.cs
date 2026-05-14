using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// Tests for iterative PEP training in <see cref="PepAnalysisEngine"/>: the loop that, after the
    /// first cross-fit fit, re-selects positive training examples from PEP-derived q-values and
    /// retrains until the positive pool stabilizes or a max iteration count is reached.
    ///
    /// These tests pin five contracts: (1) the iteration count is bounded by a configurable max,
    /// (2) the positive pool size is logged per iteration, (3) convergence stops the loop early,
    /// (4) with iterative training off the engine runs exactly one pass and never relabels (the
    /// regression safety net), and (5) cross-fit isolation - fold k is relabeled only from the
    /// model trained on folds != k, never from a single global PEP.
    /// </summary>
    [TestFixture]
    public static class PepIterativeTrainingTest
    {
        private const string MzmlName = "TaGe_SA_HeLa_04_subset_longestSeq.mzML";
        private static string OutputFolder => Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\");

        private static List<(string, CommonParameters)> _fsp;
        private static MassSpectrometry.MsDataFile _msDataFile;
        private static List<Protein> _proteinList;
        private static CommonParameters _commonParameters;

        /// <summary>
        /// Loads the spectra file and protein database once. The actual search is re-run per test
        /// (see <see cref="BuildFreshPsms"/>) because the PEP engine mutates the PSMs it scores.
        /// </summary>
        [OneTimeSetUp]
        public static void SetUp()
        {
            _commonParameters = new CommonParameters(digestionParams: new DigestionParams());
            _fsp = new List<(string, CommonParameters)> { (MzmlName, _commonParameters) };

            string mzmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\" + MzmlName);
            _msDataFile = new MyFileManager(true).LoadFile(mzmlPath, _commonParameters);

            _proteinList = ProteinDbLoader.LoadProteinFasta(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\hela_snip_for_unitTest.fasta"),
                true, DecoyType.Reverse, false, out _,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex,
                ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex, -1);
        }

        /// <summary>
        /// Runs a fresh classic search + FDR analysis and returns the resulting PSMs with raw-score
        /// q-values populated (the input state the PEP engine expects). Each call produces independent
        /// PSM objects so tests do not contaminate one another through the engine's mutations.
        /// </summary>
        private static List<SpectralMatch> BuildFreshPsms()
        {
            var ms2Scans = MetaMorpheusTask.GetMs2Scans(_msDataFile, @"TestData\" + MzmlName, _commonParameters)
                .OrderBy(b => b.PrecursorMass).ToArray();

            SpectralMatch[] allPsmsArray = new PeptideSpectralMatch[ms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, ms2Scans, new List<Modification>(), new List<Modification>(),
                null, null, null, _proteinList, new SinglePpmAroundZeroSearchMode(5), _commonParameters, _fsp,
                null, new List<string>(), false).Run();

            var psms = allPsmsArray.Where(p => p != null).ToList();
            // populates FdrInfo.QValue (the raw-score q-value used to seed iteration-0 labels)
            new FdrAnalysisEngine(psms, 1, _commonParameters, _fsp, new List<string>()).Run();
            return psms;
        }

        private static PepAnalysisEngine MakeEngine(List<SpectralMatch> psms, bool iterativeTraining, int maxTrainingIterations = 3)
        {
            return new PepAnalysisEngine(psms, "standard", _fsp, OutputFolder,
                iterativeTraining: iterativeTraining, maxTrainingIterations: maxTrainingIterations);
        }

        // ---- 1. Iteration count is bounded by a configurable max -----------------------------------

        [Test]
        public static void MaxTrainingIterations_IsClampedToAtLeastOne()
        {
            // The max is configurable via the constructor; values below 1 are nonsensical and clamped.
            var psms = BuildFreshPsms();
            Assert.That(MakeEngine(psms, iterativeTraining: true, maxTrainingIterations: 5).MaxTrainingIterations, Is.EqualTo(5));
            Assert.That(MakeEngine(psms, iterativeTraining: true, maxTrainingIterations: 1).MaxTrainingIterations, Is.EqualTo(1));
            Assert.That(MakeEngine(psms, iterativeTraining: true, maxTrainingIterations: 0).MaxTrainingIterations, Is.EqualTo(1));
            Assert.That(MakeEngine(psms, iterativeTraining: true, maxTrainingIterations: -5).MaxTrainingIterations, Is.EqualTo(1));
        }

        [Test]
        public static void IterativeTraining_RunsNoMoreThanMaxIterations()
        {
            var psms = BuildFreshPsms();
            var engine = MakeEngine(psms, iterativeTraining: true, maxTrainingIterations: 2);
            string metrics = engine.ComputePEPValuesForAllPSMs();

            Assert.That(metrics, Does.Not.Contain("failed"), "test dataset was too small to train the PEP model");
            Assert.That(engine.PositivePoolCountPerIteration.Count, Is.GreaterThanOrEqualTo(1));
            Assert.That(engine.PositivePoolCountPerIteration.Count, Is.LessThanOrEqualTo(2),
                "iteration count must never exceed MaxTrainingIterations");
        }

        // ---- 2. Positive pool size is logged per iteration ----------------------------------------

        [Test]
        public static void PositivePoolCount_IsRecordedOncePerIteration()
        {
            var psms = BuildFreshPsms();
            var engine = MakeEngine(psms, iterativeTraining: true, maxTrainingIterations: 3);
            string metrics = engine.ComputePEPValuesForAllPSMs();

            Assert.That(metrics, Does.Not.Contain("failed"), "test dataset was too small to train the PEP model");
            // every iteration that actually trained contributes exactly one positive-pool count
            Assert.That(engine.PositivePoolCountPerIteration, Is.Not.Empty);
            Assert.That(engine.PositivePoolCountPerIteration, Has.All.GreaterThan(0));
            // and the per-iteration counts are surfaced in the human-readable metrics output
            Assert.That(metrics, Does.Contain("Positive Training Pool Per Iteration"));
            Assert.That(metrics, Does.Contain("Training Iterations Run:  " + engine.PositivePoolCountPerIteration.Count));
        }

        // ---- 3. Convergence stops the loop early --------------------------------------------------

        [Test]
        public static void HasConverged_TriggersBelowOnePercentChange()
        {
            // < 1% change in the positive pool counts as converged; >= 1% does not.
            Assert.That(PepAnalysisEngine.HasConverged(10000, 10000), Is.True, "no change is converged");
            Assert.That(PepAnalysisEngine.HasConverged(10000, 10099), Is.True, "0.99% change is converged");
            Assert.That(PepAnalysisEngine.HasConverged(10000, 9950), Is.True, "a small decrease is converged");
            Assert.That(PepAnalysisEngine.HasConverged(10000, 10100), Is.False, "exactly 1% change is not converged");
            Assert.That(PepAnalysisEngine.HasConverged(10000, 10500), Is.False, "5% change is not converged");
            // guards against divide-by-zero / degenerate previous counts
            Assert.That(PepAnalysisEngine.HasConverged(0, 0), Is.True);
        }

        // ---- 4. Iterative training off => single pass, no relabeling (regression safety net) ------

        [Test]
        public static void IterativeTrainingOff_RunsExactlyOnePassAndNeverRelabels()
        {
            var psms = BuildFreshPsms();
            var engine = MakeEngine(psms, iterativeTraining: false);
            string metrics = engine.ComputePEPValuesForAllPSMs();

            Assert.That(metrics, Does.Not.Contain("failed"), "test dataset was too small to train the PEP model");
            Assert.That(engine.IterativeTraining, Is.False);
            Assert.That(engine.PositivePoolCountPerIteration.Count, Is.EqualTo(1),
                "with iterative training off the loop must run exactly one pass");
            Assert.That(engine.PerFoldRelabelPepQValues, Is.Null,
                "with iterative training off no relabeling map is ever built");
        }

        [Test]
        public static void IterativeTrainingOn_WithMaxOneIteration_BehavesLikeSinglePass()
        {
            // maxTrainingIterations == 1 means the very first pass is also the final pass: it labels
            // from raw q-values and never relabels, exactly as the single-pass pipeline does.
            var psms = BuildFreshPsms();
            var engine = MakeEngine(psms, iterativeTraining: true, maxTrainingIterations: 1);
            string metrics = engine.ComputePEPValuesForAllPSMs();

            Assert.That(metrics, Does.Not.Contain("failed"), "test dataset was too small to train the PEP model");
            Assert.That(engine.PositivePoolCountPerIteration.Count, Is.EqualTo(1));
            Assert.That(engine.PerFoldRelabelPepQValues, Is.Null);
        }

        // ---- 5. Cross-fit isolation ---------------------------------------------------------------

        [Test]
        public static void Relabeling_IsPerFoldAndDisjoint_NotGlobal()
        {
            // The most important contract: on iterations >= 2, fold k's positive examples are
            // re-selected from a PEP q-value map built ONLY from fold k's held-out predictions -
            // i.e. from the model trained on folds != k. A single global PEP map would instead
            // contain every PSM. This test proves the maps partition the PSMs exactly by fold.
            var psms = BuildFreshPsms();
            var engine = MakeEngine(psms, iterativeTraining: true, maxTrainingIterations: 2);
            string metrics = engine.ComputePEPValuesForAllPSMs();

            Assert.That(metrics, Does.Not.Contain("failed"), "test dataset was too small to train the PEP model");
            Assert.That(engine.PerFoldRelabelPepQValues, Is.Not.Null,
                "a two-iteration run must build the relabel maps for iteration 2");
            Assert.That(engine.PerFoldRelabelPepQValues.Count, Is.EqualTo(4), "one map per cross-fit fold");

            // Reconstruct the engine's folds from its (deterministic) public grouping + partitioning.
            var groups = engine.UsePeptideLevelQValueForTraining
                ? SpectralMatchGroup.GroupByBaseSequence(engine.AllPsms)
                : SpectralMatchGroup.GroupByIndividualPsm(engine.AllPsms);
            var foldIndices = PepAnalysisEngine.GetPeptideGroupIndices(groups, 4);

            var seenAcrossFolds = new HashSet<SpectralMatch>();
            for (int k = 0; k < 4; k++)
            {
                var expectedFoldMatches = foldIndices[k]
                    .SelectMany(idx => groups[idx].GetBestMatches())
                    .Where(psm => psm != null)
                    .ToHashSet();

                var relabelKeys = engine.PerFoldRelabelPepQValues[k].Keys.ToHashSet();

                // fold k's relabel map is keyed by exactly fold k's matches - nothing global
                Assert.That(relabelKeys, Is.EquivalentTo(expectedFoldMatches),
                    $"fold {k}'s relabel map must be keyed only by fold {k}'s PSMs");

                // and the per-fold maps never overlap - each PSM is relabeled by its own fold's model
                foreach (var match in relabelKeys)
                {
                    Assert.That(seenAcrossFolds.Add(match), Is.True,
                        "a PSM appeared in more than one fold's relabel map");
                }

                // q-values are real probabilities
                Assert.That(engine.PerFoldRelabelPepQValues[k].Values, Has.All.InRange(0.0, 1.0));
            }
        }
    }
}
