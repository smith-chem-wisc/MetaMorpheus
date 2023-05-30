using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.SpectralRecovery;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.DirectoryServices;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using Chemistry;
using Easy.Common.Extensions;
using FlashLFQ;
using TaskLayer;


namespace Test
{
    [TestFixture]
    public static class
    SpectralRecoveryTest
    {
        private static MyTaskResults searchTaskResults;
        private static List<PsmFromTsv> tsvPsms;
        private static List<PeptideSpectralMatch> psms;
        private static List<Protein> proteinList;
        private static MyFileManager myFileManager;
        private static List<string> rawSlices;
        private static List<DbForTask> databaseList;
        private static string outputFolder;
        private static string narrowWindowOutputFolder;
        private static Dictionary<string, int[]> numSpectraPerFile;

        [OneTimeSetUp]
        public static void SpectralRecoveryTestSetup()
        {
            // This block of code converts from PsmFromTsv to PeptideSpectralMatch objects
            string psmtsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"SpectralRecoveryTest\MSMSids.psmtsv");
            tsvPsms = PsmTsvReader.ReadTsv(psmtsvPath, out var warnings);
            psms = new List<PeptideSpectralMatch>();
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
                PeptideSpectralMatch psm = new PeptideSpectralMatch(pwsm, 0, readPsm.Score, readPsm.Ms2ScanNumber, ms2Scan,
                    new CommonParameters(), readPsm.MatchedIons);
                psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
                if (readPsm.Ms2ScanNumber == 206 && readPsm.BaseSeq.Equals("HADIVTTTTHK")) psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0.0046, 0); // Necessary for to be implemented "original pep" test
                psms.Add(psm);
                proteinList.Add(protein);
            }

            outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSpectralRecoveryOutput");
            narrowWindowOutputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"smallWindowSpectralRecoveryOutput");
            Directory.CreateDirectory(outputFolder);
            Directory.CreateDirectory(narrowWindowOutputFolder);

            numSpectraPerFile = new Dictionary<string, int[]> { { "K13_02ng_1min_frac1", new int[] { 8, 8 }
                }, { "K13_20ng_1min_frac1", new int[] { 8, 8 } } };
            rawSlices = new List<string> {
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"SpectralRecoveryTest\K13_02ng_1min_frac1.mzML"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"SpectralRecoveryTest\K13_20ng_1min_frac1.mzML") };
            databaseList = new List<DbForTask>() {new DbForTask(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"SpectralRecoveryTest\HumanFastaSlice.fasta"), false) };

            SearchTask searchTask = new SearchTask
            {
                SearchParameters = new SearchParameters()
                {
                    DoQuantification = false,
                    WriteSpectralLibrary = false,
                    MatchBetweenRuns = false,
                    DoSpectralRecovery = false,
                    WriteMzId = false
                },
                CommonParameters = new CommonParameters()
            };
            searchTaskResults = searchTask.RunTask(outputFolder, databaseList, rawSlices, "name");
            searchTaskResults = searchTask.RunTask(narrowWindowOutputFolder, databaseList, rawSlices, "name");

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
                        DoQuantification = true,
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

            postSearchTask.Parameters.OutputFolder = narrowWindowOutputFolder;
            postSearchTask.Parameters.SearchParameters.SpectralRecoveryWindowHalfWidth = 0.1;
            postSearchTask.Run();
        }

        [Test]
        public static void SpectralRecoveryPostSearchAnalysisTest()
        {
            string mbrAnalysisPath = Path.Combine(outputFolder, @"SpectralRecovery\RecoveredSpectra.psmtsv");
            string expectedHitsPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"SpectralRecoveryTest\ExpectedMBRHits.psmtsv");
            List<PsmFromTsv> mbrPsms = PsmTsvReader.ReadTsv(mbrAnalysisPath, out var warnings);
            // These PSMS were found in a search and removed from the MSMSids file. Theoretically, all peaks present in this file should be found by MbrAnalysis
            List<PsmFromTsv> expectedMbrPsms = PsmTsvReader.ReadTsv(expectedHitsPath, out warnings);

            List<PsmFromTsv> matches2ng = mbrPsms.Where(p => p.FileNameWithoutExtension == "K13_20ng_1min_frac1").ToList();
            List<PsmFromTsv> matches02ng = mbrPsms.Where(p => p.FileNameWithoutExtension == "K13_02ng_1min_frac1").ToList();
            List<string> expectedMatches = mbrPsms.Select(p => p.BaseSeq).Intersect(expectedMbrPsms.Select(p => p.BaseSeq).ToList()).ToList();

            // These assertions were written when spectral recovery depended on successful deconvolution to 
            // consider a recovered spectra. It has since changed to consider any scan where the isolation window
            // contains the precursor. This results in significantly more matches.
            Assert.That(matches2ng.Count >= 2);
            Assert.That(matches02ng.Count >= 8);
            Assert.That(expectedMatches.Count >= 3); // FlashLFQ doesn't find all 6 expected peaks, only 3. MbrAnalysis finds these three peaks

            // A spectrum is recovered for this peptide with the default spectral recovery window (1 minute).
            // The ms2 scan is collected 0.15 minutes away from the RT apex
            string testPeptideBaseSeq = "QVEPPAK";
            Assert.True(matches02ng
                .Any(p => p.BaseSeq.Equals(testPeptideBaseSeq)));
            // The narrow window test should only search for ms2 scans within 0.1 minutes of the RT apex, so the test
            // peptide should not be reported
            string mbrAnalysisShortWindowPath = Path.Combine(narrowWindowOutputFolder, @"SpectralRecovery\RecoveredSpectra.psmtsv");
            mbrPsms = PsmTsvReader.ReadTsv(mbrAnalysisShortWindowPath, out warnings);
            Assert.False(mbrPsms
                .Where(p => p.FileNameWithoutExtension.Equals("K13_02ng_1min_frac1"))
                .Any(p => p.BaseSeq.Equals(testPeptideBaseSeq)));

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
                    // Updated spectral recovery considers a scan that was missed in the old version,
                    // which is why the spectral contrast angle increases.
                    if (rowSplit[1].Equals("EGERPAR"))
                    {
                        Assert.That(contrastAngle, Is.EqualTo(0.7957).Within(0.001));
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
                        if (double.TryParse(rowSplit[9], out var contrastAngle))
                        {
                            Assert.That(contrastAngle, Is.EqualTo(0.7957).Within(0.001));
                        }
                        else
                        {
                            Assert.Fail();
                        }
                        
                        break;
                    }
                    else
                    {
                        Assert.That(rowSplit[7].Equals("MSMS"));
                        if (double.TryParse(rowSplit[10], out var contrastAngle))
                        {
                            Assert.That(contrastAngle, Is.EqualTo(0.7957).Within(0.001));
                        }
                        else
                        {
                            Assert.Fail();
                        }

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

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
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

            MiniClassicSearchEngine mcse = new MiniClassicSearchEngine(
                myMsDataFile,
                sl,
                commonParameters,
                null,
                myFile);

            foreach (PeptideSpectralMatch psm in allPsmsArray.Where(p => p != null))
            {
                PeptideWithSetModifications pwsm = psm.BestMatchingPeptides.First().Peptide;

                SpectralRecoveryPSM[] peptideSpectralMatches =
                    mcse.SearchAroundPeak(pwsm, null, allPsmsArray[5].ScanRetentionTime, allPsmsArray[5].ScanPrecursorCharge).ToArray();

                Assert.AreEqual(allPsmsArray[5].BaseSequence, peptideSpectralMatches[0].BaseSequence);
                // Same assertion as above, just used ToString method. For spectralMatchPSMs, we insert 10 additional
                // scores before the peptide info is given.
                Assert.AreEqual(allPsmsArray[5].ToString().Split('\t')[12], peptideSpectralMatches[0].ToString().Split('\t')[22]);
                Assert.That(peptideSpectralMatches[0].SpectralAngle, Is.EqualTo(allPsmsArray[5].SpectralAngle).Within(0.01));
            }

            PeptideWithSetModifications fakePeptide = new PeptideWithSetModifications(
                "ABCDEFG", null, p: new("ABCDEFG", "FAKEPROT"), digestionParams: new DigestionParams(),
                oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 7);
            Assert.Null(mcse.SearchAroundPeak(fakePeptide, peakRetentionTime: 31.2, peakCharge: 2));
        }

        [Test]
        public static void TestPrecursorIsotopeCorrelation()
        {
            ModificationMotif.TryGetMotif("G", out var modMotif);
            Dictionary<string, Modification> modDict =
                new()
                {
                    { "Carbamidomethyl on C", GlobalVariables.AllModsKnownDictionary["Carbamidomethyl on C"] },
                    { "BS on G" , new Modification(_originalId: "BS on G", _modificationType: "BS", _target: modMotif, _monoisotopicMass: 96.0875)}
                };
            PeptideWithSetModifications pwsm_Mods = new PeptideWithSetModifications(
                "ENQGDETQG[Speculative:BS on G]C[Common Fixed:Carbamidomethyl on C]PPQR", modDict, p: new Protein("ENQGDETQGCPPQR", "FakeProtein"), digestionParams: new DigestionParams(),
                oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 14);

            // this pwsm has an equivalent mass to the pwsm with mods, but a different chemical formula.
            PeptideWithSetModifications pwsm_EqualMass = new PeptideWithSetModifications(
                "ENQGDETQGQQPPQR", null, p: new Protein("ENQGDETQGQQPPQR", "FakeProtein2"), digestionParams: new DigestionParams(),
                oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 14);

            // Create the ChromatographicPeak
            SpectraFileInfo fakeFileInfo =
                new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "DoesNotExist.raw"), "A", 1, 1, 1);
            Identification id = new Identification(fakeFileInfo, pwsm_Mods.BaseSequence, pwsm_Mods.FullSequence,
                pwsm_Mods.MonoisotopicMass, ms2RetentionTimeInMinutes: 1.01, 1, new List<FlashLFQ.ProteinGroup>());
            ChromatographicPeak peak = new ChromatographicPeak(id, false, fakeFileInfo);

            // Get theoretical isotopic distribution from equal mass peptide
            ChemicalFormula equalMassFormula = new Proteomics.AminoAcidPolymer.Peptide(pwsm_EqualMass.BaseSequence).GetChemicalFormula();
            IsotopicDistribution peptideDistribution = IsotopicDistribution
                .GetDistribution(equalMassFormula, fineResolution: 0.01, minProbability: 0.005);
            double[] theoreticalMasses = peptideDistribution.Masses;
            double[] theoreticalIntensities = peptideDistribution.Intensities;
            List<FlashLFQ.IsotopicEnvelope> lfqEnvelopes = new List<FlashLFQ.IsotopicEnvelope>();

            // Create IndexedPeaks and IsotopicEnvelopes
            for (int i = 0; i < peptideDistribution.Masses.Length; i++)
            {
                lfqEnvelopes.Add(
                    new FlashLFQ.IsotopicEnvelope(
                        new IndexedMassSpectralPeak(mz: theoreticalMasses[i].ToMz(1), intensity: theoreticalIntensities[i], 0, 1.00),
                        chargeState: 1,
                        intensity: theoreticalIntensities[i]));
            }

            // Set private properties for the ChromatographicPeak
            peak.IsotopicEnvelopes = lfqEnvelopes;
            peak.SetChromatographicPeakProperties("Apex", lfqEnvelopes.MaxBy(i => i.Intensity));

            double? equalMassCorrelation = SpectralRecoveryPSM.GetIsotopeCorrelation(pwsm_EqualMass, peak);
            double? modsCorrelation = SpectralRecoveryPSM.GetIsotopeCorrelation(pwsm_Mods, peak);

            Assert.That(equalMassCorrelation, Is.EqualTo(1).Within(0.001));
            Assert.That(modsCorrelation != null && modsCorrelation >= 0.9);

            // Now, testing one mod (Carbamidomethyl on C)
            modDict.Remove("BS on G");
            PeptideWithSetModifications pwsm_OneMod = new PeptideWithSetModifications(
                "ENQGDETQGC[Common Fixed:Carbamidomethyl on C]PPQR", modDict, p: new Protein("ENQGDETQGCPPQR", "FakeProtein"), digestionParams: new DigestionParams(),
                oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 14);
            // Get theoretical isotopic distribution from equal mass peptide
            equalMassFormula = pwsm_OneMod.FullChemicalFormula;
            peptideDistribution = IsotopicDistribution
                .GetDistribution(equalMassFormula, fineResolution: 0.01, minProbability: 0.005);
            theoreticalMasses = peptideDistribution.Masses;
            theoreticalIntensities = peptideDistribution.Intensities;
            lfqEnvelopes = new List<FlashLFQ.IsotopicEnvelope>();

            // Create IndexedPeaks and IsotopicEnvelopes
            for (int i = 0; i < peptideDistribution.Masses.Length; i++)
            {
                lfqEnvelopes.Add(
                    new FlashLFQ.IsotopicEnvelope(
                        new IndexedMassSpectralPeak(mz: theoreticalMasses[i].ToMz(1), intensity: theoreticalIntensities[i], 0, 1.00),
                        chargeState: 1,
                        intensity: theoreticalIntensities[i]));
            }

            // Set private properties for the ChromatographicPeak
            peak.IsotopicEnvelopes = lfqEnvelopes;
            peak.SetChromatographicPeakProperties("Apex", lfqEnvelopes.MaxBy(i => i.Intensity));
            double? oneModCorrelation = SpectralRecoveryPSM.GetIsotopeCorrelation(pwsm_OneMod, peak);

            Assert.That(oneModCorrelation, Is.EqualTo(1).Within(0.001));
        }

        public static void SetChromatographicPeakProperties(this ChromatographicPeak peak, string propName, Object newValue)
        {
            PropertyInfo propertyInfo = typeof(ChromatographicPeak).GetProperty(propName);
            if (propertyInfo == null || propertyInfo.PropertyType != newValue.GetType()) return;
            propertyInfo.SetValue(peak, newValue);
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
                        DoQuantification = true,
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
                        DoQuantification = true,
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
            string[] psmHeader = PeptideSpectralMatch.GetTabSeparatedHeader().Trim().Split('\t');

            List<string> expectedHeader = new List<string>(psmHeader[..12]);
            expectedHeader.AddRange( new List<string>{ 
                "Peak Apex RT (min)",
                "RT Shift",
                "Deconvolutable Precursor",
                "Precursor Isotopic Envelope Score",
                "Isolation Window Center (Th)",
                "Precursor m/z - Isolation Center Distance (Th)",
                "Isolation Window Width (Th)",
                "Original Psm QValue",
                "Original Psm PEP",
                "Original Psm PEP_QValue"});
            expectedHeader.AddRange(psmHeader[12..]);

            Assert.AreEqual(String.Join('\t', expectedHeader), SpectralRecoveryPSM.TabSeparatedHeader);
        }

        [OneTimeTearDown]
        public static void SpectralRecoveryTeardown()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSpectralRecoveryOutput");
            Directory.Delete(filePath, true);
        }
    }
}