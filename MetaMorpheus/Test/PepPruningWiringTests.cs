using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.FdrAnalysis;
using EngineLayer.GlycoSearch;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// Pins which searches let PEP prune ambiguous match hypotheses. Pruning is on only where no DisambiguationEngine
    /// runs afterwards: glyco, crosslink, and nonspecific or semi-specific searches (PostSearchAnalysisTask skips its
    /// FDR pass and disambiguation for SearchType.NonSpecific). A classic SearchTask must not prune.
    /// </summary>
    [TestFixture]
    [NonParallelizable]
    public static class PepPruningWiringTests
    {
        private static string OutputFolder => Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPepPruningWiring");

        [TearDown]
        public static void TearDown()
        {
            if (Directory.Exists(OutputFolder))
            {
                Directory.Delete(OutputFolder, true);
            }
        }

        /// <summary>
        /// Runs the task and returns the pruning setting of every FdrAnalysisEngine it ran, read from the engine itself,
        /// so each call site's choice is pinned whether or not the small test data is enough to train a PEP model.
        /// </summary>
        private static List<bool> FdrEnginePruningSettings(MetaMorpheusTask task, string database, List<string> spectra)
        {
            var field = typeof(FdrAnalysisEngine).GetField("PruneAmbiguousHypotheses", BindingFlags.NonPublic | BindingFlags.Instance)!;
            var settings = new List<bool>();
            void Handler(object sender, SingleEngineFinishedEventArgs e)
            {
                if (e.MyResults is FdrAnalysisResults { MyEngine: FdrAnalysisEngine engine })
                {
                    lock (settings)
                    {
                        settings.Add((bool)field.GetValue(engine)!);
                    }
                }
            }

            MetaMorpheusEngine.FinishedSingleEngineHandler += Handler;
            try
            {
                Directory.CreateDirectory(OutputFolder);
                task.RunTask(OutputFolder, new List<DbForTask> { new DbForTask(database, false) }, spectra, "wiring");
            }
            finally
            {
                MetaMorpheusEngine.FinishedSingleEngineHandler -= Handler;
            }

            // Without an FDR run the assertions below would pass vacuously.
            Assert.That(settings, Is.Not.Empty, "no FdrAnalysisEngine ran");
            return settings;
        }

        [Test]
        [TestCase(SearchType.Classic, false)]
        [TestCase(SearchType.NonSpecific, true)]
        public static void SearchTask_PrunesOnlyForNonspecificSearch(SearchType searchType, bool expectedPruning)
        {
            var task = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    SearchType = searchType,
                    LocalFdrCategories = new List<FdrCategory> { FdrCategory.FullySpecific, FdrCategory.SemiSpecific }
                },
                CommonParameters = new CommonParameters(scoreCutoff: 4, addCompIons: true,
                    digestionParams: new DigestionParams(searchModeType: searchType == SearchType.NonSpecific ? CleavageSpecificity.Semi : CleavageSpecificity.Full,
                        fragmentationTerminus: searchType == SearchType.NonSpecific ? FragmentationTerminus.N : FragmentationTerminus.Both))
            };
            var settings = FdrEnginePruningSettings(task,
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta"),
                new List<string> { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml") });

            Assert.That(settings, Is.All.EqualTo(expectedPruning));
        }

        [Test]
        public static void CrosslinkSearch_Prunes()
        {
            var task = new XLSearchTask();
            task.XlSearchParameters.Crosslinker = GlobalVariables.Crosslinkers.ToList()[1];
            var settings = FdrEnginePruningSettings(task,
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA.fasta"),
                new List<string> { Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA_DSS_23747.mzML") });

            Assert.That(settings, Is.All.True);
        }

        [Test]
        public static void GlycoSearch_Prunes()
        {
            string dataFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\N_O_glycoWithFileSpecific");
            var task = new GlycoSearchTask
            {
                CommonParameters = new CommonParameters(dissociationType: DissociationType.HCD, ms2childScanDissociationType: DissociationType.EThcD),
                _glycoSearchParameters = new GlycoSearchParameters
                {
                    OGlycanDatabasefile = "OGlycan.gdb",
                    GlycoSearchType = GlycoSearchType.OGlycanSearch,
                    OxoniumIonFilt = true,
                    DecoyType = DecoyType.Reverse,
                    GlycoSearchTopNum = 50,
                    MaximumOGlycanAllowed = 4,
                    DoParsimony = true,
                    WritePrunedDataBase = false,
                }
            };
            var settings = FdrEnginePruningSettings(task,
                Path.Combine(dataFolder, "FourMucins_NoSigPeps_FASTA.fasta"),
                Directory.GetFiles(dataFolder).Where(p => p.Contains("mzML")).ToList());

            Assert.That(settings, Is.All.True);
        }
    }
}
