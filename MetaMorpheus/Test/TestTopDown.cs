using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using Easy.Common.Extensions;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using iText.Forms.Xfdf;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Fragmentation.Peptide;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using TaskLayer;
using UsefulProteomicsDatabases;
using Mzml = IO.MzML.Mzml;

namespace Test
{
    [TestFixture]
    public class TestTopDown
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
            Assert.That(GlobalVariables.AnalyteType == "Proteoform");

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

            SpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            bool writeSpetralLibrary = false;
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null, 
                proteinList, searchMode, CommonParameters, null, null, new List<string>(), writeSpetralLibrary).Run();

            var psm = allPsmsArray.Where(p => p != null).FirstOrDefault();
            Assert.That(psm.MatchedFragmentIons.Count == 47);
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

            SpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 1, DecoyType.Reverse, CommonParameters,
                null, 30000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();

            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, 0, CommonParameters, null, searchMode, 0, new List<string>()).Run();

            var psm = allPsmsArray.Where(p => p != null).FirstOrDefault();
            Assert.That(psm.MatchedFragmentIons.Count == 47);
        }

        [Test]
        public static void TestTopNIdentificationsPerSpectrum()
        {
            string tempDirPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "ChimeraFilteringTest");
            if (!Directory.Exists(tempDirPath))
                Directory.CreateDirectory(tempDirPath);

            var spectraPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData",
                "HighlyChimericSnip.mzML");
            var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData",
                "HistoneH3GPTMD.xml");
            var taskPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData",
                               "HighlyChimeric.toml");
            var searchTask = Toml.ReadFile<SearchTask>(taskPath, MetaMorpheusTask.tomlConfig);
            
            var maxIdsToTry = new int[] { 100, 10, 5, 1 };
            foreach (var maxIds in maxIdsToTry)
            {
                searchTask.CommonParameters.MaximumIdentificationsPerSpectrum = maxIds;
                var searchTaskDir = Path.Combine(tempDirPath, $"Task1-SearchTask{maxIds}");
                var runner = new EverythingRunnerEngine(
                    new List<(string, MetaMorpheusTask)> { ($"Task1-SearchTask{maxIds}", searchTask) },
                    new List<string> { spectraPath },
                    new List<DbForTask> { new DbForTask(dbPath, false) }, tempDirPath);
                runner.Run();


                var psms = PsmTsvReader.ReadTsv(Path.Combine(searchTaskDir, "AllPSMs.psmtsv"), out _);
                    
                    psms.GroupBy(psm => psm, new CustomComparer<PsmFromTsv>(psm => psm.Ms2ScanNumber, psm => psm.FileNameWithoutExtension))
                    .ForEach(chimeraGroup => Assert.That(chimeraGroup.Count(), Is.LessThanOrEqualTo(maxIds)));

                var proteoforms = PsmTsvReader.ReadTsv(Path.Combine(searchTaskDir, "AllProteoforms.psmtsv"), out _);
                    proteoforms.GroupBy(proteoform => proteoform, CustomComparer<PsmFromTsv>.ChimeraComparer)
                    .ForEach(chimeraGroup => Assert.That(chimeraGroup.Count(), Is.LessThanOrEqualTo(maxIds)));
            }

            Directory.Delete(tempDirPath, true);
        }
    }
}
