using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.FdrAnalysis;
using EngineLayer.GlycoSearch;
using EngineLayer.SpectrumMatch;
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
    /// Pins which searches disambiguate by PEP. The PEP engine removes no hypotheses; a DisambiguationEngine with a PEP
    /// rule does, right after each FDR run that trains PEP. Glyco, crosslink, and nonspecific or semi-specific searches
    /// use the rule the PEP engine used to apply itself. A classic SearchTask's DisambiguationEngine has no PEP rule yet.
    /// </summary>
    [TestFixture]
    [NonParallelizable]
    public static class PepDisambiguationWiringTests
    {
        private static string OutputFolder => Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPepDisambiguationWiring");

        [TearDown]
        public static void TearDown()
        {
            if (Directory.Exists(OutputFolder))
            {
                Directory.Delete(OutputFolder, true);
            }
        }

        /// <summary>
        /// One entry per FdrAnalysisEngine or DisambiguationEngine the task ran, in order: for an FDR run, whether it
        /// trains PEP (read from the engine, so the call site's choice is pinned whether or not the small test data is
        /// enough to train a model); for a disambiguation run, whether it had a PEP rule.
        /// </summary>
        private static List<(bool IsFdr, bool PepFlag)> FdrAndDisambiguationRuns(MetaMorpheusTask task, string database, List<string> spectra)
        {
            var doPep = typeof(FdrAnalysisEngine).GetField("DoPEP", BindingFlags.NonPublic | BindingFlags.Instance)!;
            var runs = new List<(bool IsFdr, bool PepFlag)>();
            void Handler(object sender, SingleEngineFinishedEventArgs e)
            {
                lock (runs)
                {
                    if (e.MyResults is FdrAnalysisResults { MyEngine: FdrAnalysisEngine engine })
                        runs.Add((true, (bool)doPep.GetValue(engine)!));
                    else if (e.MyResults is DisambiguationEngineResults disambiguation)
                        runs.Add((false, disambiguation.RemovedByPEP.HasValue));
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
            Assert.That(runs.Any(r => r.IsFdr), "no FdrAnalysisEngine ran");
            return runs;
        }

        /// <summary>
        /// Every FDR run that trains PEP is followed at once by disambiguation with a PEP rule.
        /// </summary>
        private static void AssertEveryPepRunIsDisambiguatedByPep(List<(bool IsFdr, bool PepFlag)> runs)
        {
            Assert.That(runs.Any(r => r is (true, true)), "no FDR run trains PEP");
            for (int i = 0; i < runs.Count; i++)
            {
                if (runs[i] is (true, true))
                {
                    Assert.That(i + 1 < runs.Count && runs[i + 1] == (false, true),
                        $"FDR run {i} trains PEP but is not followed by disambiguation by PEP");
                }
            }
        }

        [Test]
        public static void ClassicSearch_HasNoPepRuleYet()
        {
            var runs = FdrAndDisambiguationRuns(SearchTaskOf(SearchType.Classic),
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta"),
                new List<string> { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml") });

            Assert.That(runs.Where(r => !r.IsFdr), Is.Not.Empty, "no DisambiguationEngine ran");
            Assert.That(runs.Where(r => !r.IsFdr).Select(r => r.PepFlag), Is.All.False);
        }

        [Test]
        public static void NonspecificSearch_DisambiguatesByPep()
        {
            var runs = FdrAndDisambiguationRuns(SearchTaskOf(SearchType.NonSpecific),
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta"),
                new List<string> { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml") });

            AssertEveryPepRunIsDisambiguatedByPep(runs);
        }

        [Test]
        public static void CrosslinkSearch_DisambiguatesByPep()
        {
            var task = new XLSearchTask();
            task.XlSearchParameters.Crosslinker = GlobalVariables.Crosslinkers.ToList()[1];
            var runs = FdrAndDisambiguationRuns(task,
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA.fasta"),
                new List<string> { Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA_DSS_23747.mzML") });

            // crosslinks, singles, and the single / loop / deadend passes: all train PEP
            Assert.That(runs.Count(r => r.IsFdr), Is.EqualTo(5));
            AssertEveryPepRunIsDisambiguatedByPep(runs);
        }

        [Test]
        public static void GlycoSearch_DisambiguatesByPep()
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
            var runs = FdrAndDisambiguationRuns(task,
                Path.Combine(dataFolder, "FourMucins_NoSigPeps_FASTA.fasta"),
                Directory.GetFiles(dataFolder).Where(p => p.Contains("mzML")).ToList());

            AssertEveryPepRunIsDisambiguatedByPep(runs);
        }

        private static SearchTask SearchTaskOf(SearchType searchType) => new SearchTask
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
    }
}
