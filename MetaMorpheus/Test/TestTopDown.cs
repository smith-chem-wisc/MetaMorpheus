using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using Easy.Common.Extensions;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using NUnit.Framework;
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
                new Protein(
                    "MPKVYSYQEVAEHNGPENFWIIIDDKVYDVSQFKDEHPGGDEIIMDLGGQDATESFVDIGHSDEALRLLKGLYIGDVDKTSERVSVEKVSTSENQSKGSGTLVVILAILMLGVAYYLLNE",
                    "P40312")
            };

            var myMsDataFile = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TopDownTestData\slicedTDYeast.mzML"));

            var searchMode = new SinglePpmAroundZeroSearchMode(5);

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters())
                .OrderBy(b => b.PrecursorMass).ToArray();

            SpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            bool writeSpetralLibrary = false;
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null,
                null, null,
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
                new Protein(
                    "MPKVYSYQEVAEHNGPENFWIIIDDKVYDVSQFKDEHPGGDEIIMDLGGQDATESFVDIGHSDEALRLLKGLYIGDVDKTSERVSVEKVSTSENQSKGSGTLVVILAILMLGVAYYLLNE",
                    "P40312")
            };

            var myMsDataFile = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TopDownTestData\slicedTDYeast.mzML"));

            var searchMode = new SinglePpmAroundZeroSearchMode(5);

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters())
                .OrderBy(b => b.PrecursorMass).ToArray();

            SpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null,
                null, 1, DecoyType.Reverse, CommonParameters,
                null, 30000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant,
                new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();

            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex,
                indexResults.FragmentIndex, 0, CommonParameters, null, searchMode, 0, new List<string>()).Run();

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
            var secondSpectraPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData",
                "TDGPTMDSearchSingleSpectra.mzML");
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
                    new List<string> { spectraPath, secondSpectraPath },
                    new List<DbForTask> { new DbForTask(dbPath, false) }, tempDirPath);
                runner.Run();

                // check bulk results
                PsmTsvReader.ReadTsv(Path.Combine(searchTaskDir, "AllPSMs.psmtsv"), out _)
                    .GroupBy(psm => psm,
                        new CustomComparer<PsmFromTsv>(psm => psm.Ms2ScanNumber, psm => psm.FileNameWithoutExtension))
                    .ForEach(chimeraGroup => Assert.That(chimeraGroup.Count(), Is.LessThanOrEqualTo(maxIds)));

                PsmTsvReader.ReadTsv(Path.Combine(searchTaskDir, "AllProteoforms.psmtsv"), out _)
                    .GroupBy(proteoform => proteoform, CustomComparer<PsmFromTsv>.ChimeraComparer)
                    .ForEach(chimeraGroup => Assert.That(chimeraGroup.Count(), Is.LessThanOrEqualTo(maxIds)));

                // check individual file results
                var files = Directory.GetFiles(Path.Combine(searchTaskDir, "Individual File Results"));
                files.Where(p => p.Contains("PSMs.psmtsv"))
                    .ForEach(psmFile => PsmTsvReader.ReadTsv(psmFile, out _)
                        .GroupBy(psm => psm, CustomComparer<PsmFromTsv>.ChimeraComparer)
                        .ForEach(chimeraGroup => Assert.That(chimeraGroup.Count(), Is.LessThanOrEqualTo(maxIds))));

                files.Where(p => p.Contains("Proteoforms.psmtsv"))
                    .ForEach(proteoformFile => PsmTsvReader.ReadTsv(proteoformFile, out _)
                        .GroupBy(proteoform => proteoform, CustomComparer<PsmFromTsv>.ChimeraComparer)
                        .ForEach(chimeraGroup => Assert.That(chimeraGroup.Count(), Is.LessThanOrEqualTo(maxIds))));
            }

            Directory.Delete(tempDirPath, true);
        }


        [Test]
        public static void TestCustomComparer_PsmFromTsv()
        {
            var resultPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData",
                "HighlyChimeric.psmtsv");
            
            var psms = PsmTsvReader.ReadTsv(resultPath, out var errors);
            Assert.That(!errors.Any());

            foreach (var chimeraGroup in psms.GroupBy(p => p, CustomComparer<PsmFromTsv>.ChimeraComparer))
            {
                Assert.That(chimeraGroup.All(p => p.Ms2ScanNumber == chimeraGroup.Key.Ms2ScanNumber));
                Assert.That(chimeraGroup.All(p => p.PrecursorScanNum == chimeraGroup.Key.PrecursorScanNum));
                Assert.That(chimeraGroup.All(p => p.FileNameWithoutExtension == chimeraGroup.Key.FileNameWithoutExtension));
            }
        }

        [Test]
        public static void TestCustomComparerAndExcessChimeraRemover_SpectralMatch()
        {
            // setup
            string tempDirPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "ChimeraFilteringTest");
            if (!Directory.Exists(tempDirPath))
                Directory.CreateDirectory(tempDirPath);

            var spectraPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "HighlyChimericSnip.mzML");
            var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "HistoneH3GPTMD.xml");
            var taskPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "HighlyChimeric.toml");

            // build all components for the engine
            var searchTask = Toml.ReadFile<SearchTask>(taskPath, MetaMorpheusTask.tomlConfig);
            var commonParams = searchTask.CommonParameters;
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();

            var searchMode = SearchTask.GetMassDiffAcceptor(commonParams.PrecursorMassTolerance,
                searchTask.SearchParameters.MassDiffAcceptorType, searchTask.SearchParameters.CustomMdac);
            MetaMorpheusTask.DetermineAnalyteType(commonParams);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(MsDataFileReader.GetDataFile(spectraPath), null,
                commonParams).OrderBy(b => b.PrecursorMass).ToArray(); 
            SpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];

            // build and run the engine
            var engine = new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications,
                fixedModifications, null, null, null,
                ProteinDbLoader.LoadProteinXML(dbPath, true, DecoyType.Reverse, GlobalVariables.AllModsKnown, false,
                    new List<string>(), out var unknownMods), searchMode, commonParams, null, null, new List<string>(),
                false);
            engine.Run();

            var psms = allPsmsArray.Where(p => p != null).OrderByDescending(p => p).ToList();
            var chimeraGroups = psms.GroupBy(p => p, CustomComparer<SpectralMatch>.SMChimeraComparer)
                .ToDictionary(p => p.Key.ScanNumber, p => p.ToList());

            // filter to only have max n per group
            var maxIdsToTry = new[] { 10, 5, 2 };
            foreach (var maxIds in maxIdsToTry)
            {
                psms.RemoveExcessChimericIdentifications(maxIds);
                foreach (var chimeraGroup in psms.GroupBy(p => p, CustomComparer<SpectralMatch>.SMChimeraComparer)
                             .Select(p => p.ToArray()))
                {
                    // ensure grouping is correct
                    Assert.That(chimeraGroup.Length <= maxIds);
                    Assert.That(chimeraGroup.All(p => p.ScanNumber == chimeraGroup.First().ScanNumber));
                    Assert.That(chimeraGroup.All(p => p.PrecursorScanNumber == chimeraGroup.First().PrecursorScanNumber));
                    Assert.That(chimeraGroup.All(p => p.FullFilePath == chimeraGroup.First().FullFilePath));

                    // ensure the top n are the correct ones
                    var expectedIdentifications =
                        chimeraGroups[chimeraGroup.First().ScanNumber].OrderByDescending(p => p).Take(maxIds).ToArray();
                    Assert.AreEqual(expectedIdentifications.Length, chimeraGroup.Count());
                    for (int i = 0; i < expectedIdentifications.Length; i++)
                    {
                        var expected = expectedIdentifications[i];
                        var actual = chimeraGroup[i];
                        Assert.That(expected.Score, Is.EqualTo(actual.Score));
                        Assert.That(expected.DeltaScore, Is.EqualTo(actual.DeltaScore));
                        Assert.That(expected.MatchedFragmentIons.Count, Is.EqualTo(actual.MatchedFragmentIons.Count));
                        Assert.That(expected.MatchedFragmentIons.Select(p => p.NeutralTheoreticalProduct).SequenceEqual(
                                                       actual.MatchedFragmentIons.Select(p => p.NeutralTheoreticalProduct)));
                    }
                }
            }
        }
    }
}
