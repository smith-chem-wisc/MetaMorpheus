using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using IO.MzML;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public class TopDownTest
    {
        [Test]
        public static void TestClassicSearchEngineTopDown()
        {
            CommonParameters CommonParameters = new CommonParameters(
                digestionParams: new DigestionParams(protease: "top-down"),
                scoreCutoff: 1,
                assumeOrphanPeaksAreZ1Fragments: false);

            MetaMorpheusTask.DetermineAnalyteType(CommonParameters);

            // test output file name (should be proteoform and not peptide)
            Assert.That(GlobalVariables.AnalyteType == GlobalVariables.Analyte.Proteoform);

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein>
            {
                new Protein("MPKVYSYQEVAEHNGPENFWIIIDDKVYDVSQFKDEHPGGDEIIMDLGGQDATESFVDIGHSDEALRLLKGLYIGDVDKTSERVSVEKVSTSENQSKGSGTLVVILAILMLGVAYYLLNE", "P40312")
            };

            var myMsDataFile = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\slicedTDYeast.mzML"));

            var searchMode = new SinglePpmAroundZeroSearchMode(5);

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchMode, CommonParameters, null, null, new List<string>()).Run();

            var psm = allPsmsArray.Where(p => p != null).FirstOrDefault();
            Assert.That(psm.MatchedFragmentIons.Count == 47);
        }

        /// <summary>
        /// TODO: MetaMorpheus currently reports ambiguity at the PrSM level, but starts tossing things when we get to the proteoform/protein level. See issue #2061
        /// Example 1: a base seqeunce is needed for parsimony, but an ambiguous sequence means the base sequence is null.
        /// Example 2: a full sequence is needed for determining which peptides/proteoforms are unique, but ambiguous localization means the full sequence is null.
        /// </summary>
        [Test]
        public static void TestAmbiguousProteoformOutput()
        {
            CommonParameters commonParameters = new CommonParameters(
               digestionParams: new DigestionParams(protease: "top-down"),
               scoreCutoff: 1,
               useProvidedPrecursorInfo: false,
               deconvolutionMaxAssumedChargeState: 60,
               trimMsMsPeaks: false,
               listOfModsVariable: new List<(string, string)> { ("Common Variable", "Oxidation on M"), ("Common Biological", "Acetylation on K"), ("Common Biological", "Trimethylation on K") },
               listOfModsFixed: new List<(string, string)> { ("Common Fixed", "Carbamidomethyl on C") }
                    );

            SearchParameters searchParameters = new SearchParameters
            {
                DoQuantification = false
            };


            SearchTask searchTask = new SearchTask
            {
                CommonParameters = commonParameters,
                SearchParameters = searchParameters
            };

            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("task1", searchTask) };
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestProteoformAmbiguity");
            string mzmlName = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData/ProteoformAmbiguity.mzML");
            string fastaName = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData/ProteoformAmbiguity.fasta");
            // RUN!
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(fastaName, false) }, outputFolder);
            engine.Run();

            //There are 8 PrSMs, each with a unique proteoform and protein
            //check that all 8 PrSMs are reported, all 8 unique proteoforms, and all 8 proteins
            string[] prsmLines = File.ReadAllLines(Path.Combine(outputFolder, "task1/AllPrSMs.psmtsv"));
            Assert.AreEqual(prsmLines.Length, 9); //8 + header
            string[] proteoformLines = File.ReadAllLines(Path.Combine(outputFolder, "task1/AllProteoforms.psmtsv"));
            Assert.AreEqual(proteoformLines.Length, 4); //3 + header, five of the PrSMs have ambiguous full sequences, which are needed to determine individuality
            string[] proteinLines = File.ReadAllLines(Path.Combine(outputFolder, "task1/AllProteinGroups.tsv"));
            Assert.AreEqual(proteinLines.Length, 7); //6 + header, two of the PrSMs have ambiguous base sequences, which prevents their use in parsimony

            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestModernSearchEngineTopDown()
        {
            CommonParameters CommonParameters = new CommonParameters(
                digestionParams: new DigestionParams(protease: "top-down"),
                scoreCutoff: 1,
                assumeOrphanPeaksAreZ1Fragments: false);

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein>
            {
                new Protein("MPKVYSYQEVAEHNGPENFWIIIDDKVYDVSQFKDEHPGGDEIIMDLGGQDATESFVDIGHSDEALRLLKGLYIGDVDKTSERVSVEKVSTSENQSKGSGTLVVILAILMLGVAYYLLNE", "P40312")
            };

            var myMsDataFile = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\slicedTDYeast.mzML"));

            var searchMode = new SinglePpmAroundZeroSearchMode(5);

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 1, DecoyType.Reverse, CommonParameters,
                null, 30000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();

            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, 0, CommonParameters, null, searchMode, 0, new List<string>()).Run();

            var psm = allPsmsArray.Where(p => p != null).FirstOrDefault();
            Assert.That(psm.MatchedFragmentIons.Count == 47);
        }
    }
}