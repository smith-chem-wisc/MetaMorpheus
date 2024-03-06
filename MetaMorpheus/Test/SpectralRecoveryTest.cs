﻿using EngineLayer;
using EngineLayer.ClassicSearch;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.DirectoryServices;
using System.IO;
using System.Linq;
using System.Text;
using Easy.Common.Extensions;
using Omics.Modifications;
using TaskLayer;
using TaskLayer.MbrAnalysis;
using Omics;

namespace Test
{
    [TestFixture]
    public class
    SpectralRecoveryTest
    {
        private static MyTaskResults searchTaskResults;
        private static List<PsmFromTsv> tsvPsms;
        private static List<SpectralMatch> psms;
        private static List<Protein> proteinList;
        private static MyFileManager myFileManager;
        private static List<string> rawSlices;
        private static List<DbForTask> databaseList;
        private static string outputFolder;
        private static Dictionary<string, int[]> numSpectraPerFile;

        [OneTimeSetUp]
        public void SpectralRecoveryTestSetup()
        {
            // This block of code converts from PsmFromTsv to SpectralMatch objects
            string psmtsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"SpectralRecoveryTest\MSMSids.psmtsv");
            tsvPsms = PsmTsvReader.ReadTsv(psmtsvPath, out var warnings);
            psms = new List<SpectralMatch>();
            proteinList = new List<Protein>();
            myFileManager = new MyFileManager(true);

            foreach (PsmFromTsv readPsm in tsvPsms)
            {
                string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                    "TestData", "SpectralRecoveryTest", readPsm.FileNameWithoutExtension + ".mzML");
                MsDataScan scan = myFileManager.LoadFile(filePath, new CommonParameters()).GetOneBasedScan(readPsm.Ms2ScanNumber);
                Ms2ScanWithSpecificMass ms2Scan = new Ms2ScanWithSpecificMass(scan, readPsm.PrecursorMz, readPsm.PrecursorCharge,
                    filePath, new CommonParameters());
                Protein protein = new(readPsm.BaseSeq, readPsm.ProteinAccession, readPsm.OrganismName,
                    isDecoy: readPsm.DecoyContamTarget == "D",
                    isContaminant: readPsm.DecoyContamTarget == "C");
                string[] startAndEndResidues = readPsm.StartAndEndResiduesInProtein.Split(" ");
                int startResidue = Int32.Parse(startAndEndResidues[0].Trim('['));
                int endResidue = Int32.Parse(startAndEndResidues[2].Trim(']'));

                PeptideWithSetModifications pwsm = new PeptideWithSetModifications(
                    readPsm.FullSequence, null, p: protein, digestionParams: new DigestionParams(),
                    oneBasedStartResidueInProtein: startResidue, oneBasedEndResidueInProtein: endResidue);
                SpectralMatch psm = new PeptideSpectralMatch(pwsm, 0, readPsm.Score, readPsm.Ms2ScanNumber, ms2Scan,
                    new CommonParameters(), readPsm.MatchedIons);
                psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
                if (readPsm.Ms2ScanNumber == 206 && readPsm.BaseSeq.Equals("HADIVTTTTHK")) psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0.0046, 0); // Necessary for to be implemented "original pep" test
                psms.Add(psm);
                proteinList.Add(protein);
            }

            Directory.CreateDirectory(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSpectralRecoveryOutput"));

            numSpectraPerFile = new Dictionary<string, int[]> { { "K13_02ng_1min_frac1", new int[] { 8, 8 }
                }, { "K13_20ng_1min_frac1", new int[] { 8, 8 } } };
            rawSlices = new List<string> {
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"SpectralRecoveryTest\K13_02ng_1min_frac1.mzML"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"SpectralRecoveryTest\K13_20ng_1min_frac1.mzML") };
            databaseList = new List<DbForTask>() {new DbForTask(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"SpectralRecoveryTest\HumanFastaSlice.fasta"), false) };
            outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSpectralRecoveryOutput");

            SearchTask searchTask = new SearchTask
            {
                SearchParameters = new SearchParameters()
                {
                    DoLabelFreeQuantification = false,
                    WriteSpectralLibrary = false,
                    MatchBetweenRuns = false,
                    DoSpectralRecovery = false,
                    WriteMzId = false
                },
                CommonParameters = new CommonParameters()
            };
            searchTaskResults = searchTask.RunTask(outputFolder, databaseList, rawSlices, "name");


            PostSearchAnalysisTask postSearchTask = new PostSearchAnalysisTask()
            {
                Parameters = new PostSearchAnalysisParameters()
                {
                    ProteinList = proteinList,
                    AllPsms = psms,
                    CurrentRawFileList = rawSlices,
                    DatabaseFilenameList = databaseList,
                    OutputFolder = outputFolder,
                    NumMs2SpectraPerFile = numSpectraPerFile,
                    ListOfDigestionParams = new HashSet<DigestionParams> { new DigestionParams(generateUnlabeledProteinsForSilac: false) },
                    SearchTaskResults = searchTaskResults,
                    MyFileManager = myFileManager,
                    IndividualResultsOutputFolder = Path.Combine(outputFolder, "individual"),
                    SearchParameters = new SearchParameters()
                    {
                        DoLabelFreeQuantification = true,
                        WriteSpectralLibrary = true,
                        MatchBetweenRuns = true,
                        DoSpectralRecovery = true,
                        WriteMzId = false,
                        WriteDecoys = false,
                        WriteContaminants = false,
                        QuantifyPpmTol = 25
                    }
                },
                CommonParameters = new CommonParameters(dissociationType: DissociationType.Autodetect),
                FileSpecificParameters = new List<(string FileName, CommonParameters Parameters)> {
                    (rawSlices[0], new CommonParameters()),
                    (rawSlices[1], new CommonParameters())
                }
            };

            postSearchTask.Run();

        }

        [Test]
        public static void SpectralRecoveryPostSearchAnalysisTest()
        {

            List<string> warnings;
            string mbrAnalysisPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSpectralRecoveryOutput\SpectralRecovery\RecoveredSpectra.psmtsv");
            string expectedHitsPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"SpectralRecoveryTest\ExpectedMBRHits.psmtsv");
            List<PsmFromTsv> mbrPsms = PsmTsvReader.ReadTsv(mbrAnalysisPath, out warnings);
            // These PSMS were found in a search and removed from the MSMSids file. Theoretically, all peaks present in this file should be found by MbrAnalysis
            List<PsmFromTsv> expectedMbrPsms = PsmTsvReader.ReadTsv(expectedHitsPath, out warnings);

            List<PsmFromTsv> matches2ng = mbrPsms.Where(p => p.FileNameWithoutExtension == "K13_20ng_1min_frac1").ToList();
            List<PsmFromTsv> matches02ng = mbrPsms.Where(p => p.FileNameWithoutExtension == "K13_02ng_1min_frac1").ToList();
            List<string> expectedMatches = mbrPsms.Select(p => p.BaseSeq).Intersect(expectedMbrPsms.Select(p => p.BaseSeq).ToList()).ToList();

            Assert.That(matches2ng.Count >= 2);
            Assert.That(matches02ng.Count >= 8);
            Assert.That(expectedMatches.Count >= 3); // FlashLFQ doesn't find all 6 expected peaks, only 3. MbrAnalysis finds these three peaks

            //TODO: Add test for recovering fdrInfo from original. Currently, PsmTsvReader doesn't support the new columns, so it's hard to test
        }

        [Test]
        public static void TestSpectralRecoveryResultsOutput()
        {
            // Test that AllQuantifiedPeaks was written correctly
            string referenceDataPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestSpectralRecoveryOutput\AllQuantifiedPeaks.tsv");
            string[] peaksResults = File.ReadAllLines(referenceDataPath).Skip(1).ToArray();

            foreach (string row in peaksResults)
            {
                string[] rowSplit = row.Split('\t');
                if (rowSplit[15].Equals("MBR") && double.TryParse(rowSplit[16], out var contrastAngle))
                {
                    if (rowSplit[1].Equals("EGERPAR"))
                    {
                        Assert.That(contrastAngle, Is.EqualTo(0.6567).Within(0.001));
                        break;
                    }
                }
            }

            referenceDataPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestSpectralRecoveryOutput\AllQuantifiedPeptides.tsv");

            string[] peptideResults = File.ReadAllLines(referenceDataPath).Skip(1).ToArray();

            foreach (string row in peptideResults)
            {
                string[] rowSplit = row.Split('\t');
                if (rowSplit[0].Equals("EGERPAR"))
                {
                    if (rowSplit[7].Equals("MBR"))
                    {
                        Assert.That(rowSplit[8].Equals("MSMS"));
                        Assert.That(double.TryParse(rowSplit[9], out var contrastAngle) &&
                                    Math.Abs(contrastAngle - 0.6567) < 0.001);
                        break;
                    }
                    else
                    {
                        Assert.That(rowSplit[7].Equals("MSMS"));
                        Assert.That(double.TryParse(rowSplit[10], out var contrastAngle) &&
                                    Math.Abs(contrastAngle - 0.6567) < 0.001);
                        break;
                    }
                }
            }
        }

        [Test]
        public static void MiniClassicSearchEngineTest()
        {
            var testDir = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch");
            string myFile = Path.Combine(testDir, @"slicedMouse.raw");
            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters commonParameters = new CommonParameters(maxThreadsToUsePerFile: 1, scoreCutoff: 1);
            MsDataFile myMsDataFile = myFileManager.LoadFile(myFile, commonParameters);

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein> { new Protein("VIHDNFGIVEGLMTTVHAITATQK", "P16858") };

            string targetSpectralLibrary = Path.Combine(testDir, @"P16858_target.msp");
            string decoySpectralLibrary = Path.Combine(testDir, @"P16858_decoy.msp");

            List<string> specLibs = new List<string> { targetSpectralLibrary, decoySpectralLibrary };

            SpectralLibrary sl = new(specLibs);

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            Ms2ScanWithSpecificMass[] listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            SpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            bool writeSpectralLibrary = false;
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchModes, commonParameters, null, sl, new List<string>(), writeSpectralLibrary).Run();

            // Single search mode
            Assert.AreEqual(7, allPsmsArray.Length);
            Assert.IsTrue(allPsmsArray[5].Score > 38);
            Assert.AreEqual("VIHDNFGIVEGLMTTVHAITATQK", allPsmsArray[5].BaseSequence);
            Assert.IsTrue(!allPsmsArray[5].IsDecoy);

            SpectralLibrarySearchFunction.CalculateSpectralAngles(sl, allPsmsArray, listOfSortedms2Scans, commonParameters);
            Assert.That(allPsmsArray[5].SpectralAngle, Is.EqualTo(0.82).Within(0.01));

            foreach (SpectralMatch psm in allPsmsArray.Where(p => p != null))
            {
                IBioPolymerWithSetMods pwsm = psm.BestMatchingBioPolymersWithSetMods.First().Peptide;

                MiniClassicSearchEngine mcse = new MiniClassicSearchEngine(
                    listOfSortedms2Scans.OrderBy(p => p.RetentionTime).ToArray(),
                    searchModes,
                    commonParameters,
                    sl,
                    null);

                SpectralMatch[] peptideSpectralMatches =
                    mcse.SearchAroundPeak(pwsm, allPsmsArray[5].ScanRetentionTime).ToArray();

                Assert.AreEqual(allPsmsArray[5].BaseSequence, peptideSpectralMatches[0].BaseSequence);
                Assert.That(peptideSpectralMatches[0].SpectralAngle, Is.EqualTo(allPsmsArray[5].SpectralAngle).Within(0.01));
            }
        }

        [Test]
        public static void SpectralWriterTest()
        {

            PostSearchAnalysisTask postSearchTask = new PostSearchAnalysisTask()
            {
                Parameters = new PostSearchAnalysisParameters()
                {
                    ProteinList = proteinList,
                    AllPsms = psms,
                    CurrentRawFileList = rawSlices,
                    DatabaseFilenameList = databaseList,
                    OutputFolder = outputFolder,
                    NumMs2SpectraPerFile = numSpectraPerFile,
                    ListOfDigestionParams = new HashSet<DigestionParams> { new DigestionParams(generateUnlabeledProteinsForSilac: false) },
                    SearchTaskResults = searchTaskResults,
                    MyFileManager = myFileManager,
                    IndividualResultsOutputFolder = Path.Combine(outputFolder, "Individual File Results"),
                    SearchParameters = new SearchParameters()
                    {
                        DoLabelFreeQuantification = true,
                        WriteSpectralLibrary = true,
                        MatchBetweenRuns = true,
                        DoSpectralRecovery = true,
                        WriteMzId = false,
                        WriteDecoys = false,
                        WriteContaminants = false,
                        QuantifyPpmTol = 25
                    }
                },
                CommonParameters = new CommonParameters(dissociationType: DissociationType.Autodetect),
                FileSpecificParameters = new List<(string FileName, CommonParameters Parameters)> {
                    (rawSlices[0], new CommonParameters()),
                    (rawSlices[1], new CommonParameters())
                }
            };

            postSearchTask.Run();

            var path = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSpectralRecoveryOutput");
            var list = Directory.GetFiles(path, "*.*", SearchOption.AllDirectories);
            string matchingvalue = list.Where(p => p.Contains("SpectralLibrary")).First().ToString();
            var testLibraryWithoutDecoy = new SpectralLibrary(new List<string> { Path.Combine(path, matchingvalue) });

            Assert.That(testLibraryWithoutDecoy.TryGetSpectrum("IAGQVAAANK", 2, out var spectrum));
            Assert.That(testLibraryWithoutDecoy.TryGetSpectrum("HEVSASTQSTPASSR", 3, out spectrum));

            testLibraryWithoutDecoy.CloseConnections();

            // new task with less than 100 psms.
            postSearchTask = new PostSearchAnalysisTask()
            {
                Parameters = new PostSearchAnalysisParameters()
                {
                    ProteinList = proteinList,
                    AllPsms = psms.GetRange(0, 50),
                    CurrentRawFileList = rawSlices,
                    DatabaseFilenameList = databaseList,
                    OutputFolder = outputFolder,
                    NumMs2SpectraPerFile = numSpectraPerFile,
                    ListOfDigestionParams = new HashSet<DigestionParams> { new DigestionParams(generateUnlabeledProteinsForSilac: false) },
                    SearchTaskResults = searchTaskResults,
                    MyFileManager = myFileManager,
                    IndividualResultsOutputFolder = Path.Combine(outputFolder, "Individual File Results"),
                    SearchParameters = new SearchParameters()
                    {
                        DoLabelFreeQuantification = true,
                        WriteSpectralLibrary = true,
                        MatchBetweenRuns = true,
                        DoSpectralRecovery = true,
                        WriteMzId = false,
                        WriteDecoys = false,
                        WriteContaminants = false,
                        QuantifyPpmTol = 25
                    }
                },
                CommonParameters = new CommonParameters(dissociationType: DissociationType.Autodetect),
                FileSpecificParameters = new List<(string FileName, CommonParameters Parameters)> {
                    (rawSlices[0], new CommonParameters()),
                    (rawSlices[1], new CommonParameters())
                }
            };

            postSearchTask.Run();

            testLibraryWithoutDecoy = new SpectralLibrary(new List<string> { Path.Combine(path, matchingvalue) });

            Assert.That(testLibraryWithoutDecoy.TryGetSpectrum("EESGKPGAHVTVK", 2, out spectrum));

            // Test spectral library update
            postSearchTask.Parameters.SearchParameters.UpdateSpectralLibrary = true;
            postSearchTask.Parameters.SpectralLibrary = testLibraryWithoutDecoy;
            postSearchTask.Run();

            var libraryList = Directory.GetFiles(path, "*.*", SearchOption.AllDirectories);
            string updateLibraryPath = libraryList.First(p => p.Contains("SpectralLibrary") && !p.Contains(matchingvalue)).ToString();
            var updatedLibraryWithoutDecoy = new SpectralLibrary(new List<string> { Path.Combine(path, updateLibraryPath) });
            Assert.That(updatedLibraryWithoutDecoy.TryGetSpectrum("EESGKPGAHVTVK", 2, out spectrum));

            testLibraryWithoutDecoy.CloseConnections(); 
            updatedLibraryWithoutDecoy.CloseConnections();
        }

        [Test]
        public static void SpectralRecoveryHeaderTest()
        {
            string psmHeader = SpectralMatch.GetTabSeparatedHeader().Trim();
            StringBuilder sb = new();
            sb.Append(psmHeader);
            sb.Append('\t');
            sb.Append("Initial Search Q-Value");
            sb.Append('\t');
            sb.Append("Initial Search PEP");
            sb.Append('\t');
            sb.Append("Initial Search PEP Q-Value");

            Assert.AreEqual(sb.ToString(), SpectralRecoveryPSM.TabSeparatedHeader);

        }

        [OneTimeTearDown]
        public static void SpectralRecoveryTeardown()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSpectralRecoveryOutput");
            Directory.Delete(filePath, true);
        }
    }
}