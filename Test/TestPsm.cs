using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.Localization;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class TestPsm
    {
        [Test]
        public static void TestPsmHeader()
        {
            CommonParameters commonParameters = new CommonParameters();
            PeptideWithSetModifications pepWithSetMods = new Protein(
                "MQQQQQQQ",
                "accession1",
                "org",
                new List<Tuple<string, string>> { new Tuple<string, string>("geneNameType", "geneName") },
                new Dictionary<int, List<Modification>> { { 2, new List<Modification> { new Modification("mod", "mod") } } },
                name: "name",
                fullName: "fullName",
                sequenceVariations: new List<SequenceVariation> { new SequenceVariation(2, "P", "Q", "changed this sequence") })
                    .Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()).First();
            MsDataFile myMsDataFile = new TestDataFile(pepWithSetMods, "quadratic");
            MsDataScan scann = myMsDataFile.GetOneBasedScan(2);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(scann, 4, 1, null, new CommonParameters());

            var theoreticalIons = new List<Product>();
            pepWithSetMods.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theoreticalIons);
            var matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan, theoreticalIons, new CommonParameters());
            PeptideSpectralMatch psm = new PeptideSpectralMatch(pepWithSetMods, 1, 2, 3, scan, commonParameters, matchedIons);
            psm.ResolveAllAmbiguities();

            var t = psm.ToString();
            var tabsepheader = PeptideSpectralMatch.GetTabSeparatedHeader();
            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));

            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));

            Tolerance fragmentTolerance = new PpmTolerance(10);
            new LocalizationEngine(new List<PeptideSpectralMatch> { psm }, myMsDataFile, new CommonParameters(productMassTolerance: fragmentTolerance), null, new List<string>()).Run();

            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));

            psm.SetFdrValues(6, 6, 6, 6, 6, 0, 0, 0);

            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));
        }

        [Test]
        public static void TestQValueFilter()
        {
            SearchTask searchTask = new SearchTask()
            {
                CommonParameters = new CommonParameters
                (
                    qValueOutputFilter: 1
                )
            };

            SearchTask searchTask2 = new SearchTask()
            {
                CommonParameters = new CommonParameters
                (
                    qValueOutputFilter: 0
                )
            };

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPSMOutput");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");

            var engine = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("QValueTest", searchTask) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engine.Run();

            string psmFile = Path.Combine(outputFolder, @"QValueTest\AllPSMs.psmtsv");
            var lines = File.ReadAllLines(psmFile);
            Assert.That(lines.Length == 11);

            var engine2 = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("QValueTest", searchTask2) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engine2.Run();

            var lines2 = File.ReadAllLines(psmFile);
            Assert.That(lines2.Length == 7);
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestLongestFragmentIonSequence()
        {
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_HeLa_04_subset_longestSeq.mzML");
            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters CommonParameters = new CommonParameters();
            var myMsDataFile = myFileManager.LoadFile(origDataFile, CommonParameters);
            var searchModes = new SinglePpmAroundZeroSearchMode(5);
            List<Protein> proteinList = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\hela_snip_for_unitTest.fasta"), true, DecoyType.Reverse, false, out var dbErrors,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                    ProteinDbLoader.UniprotOrganismRegex, -1);
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();
            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            bool writeSpetralLibrary = false;
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null, proteinList, searchModes, 
                new CommonParameters(), null, null, new List<string>(), writeSpetralLibrary).Run();

            List<int> longestSeriesObserved = new List<int>();

            List<int> longestSeriesExpected = new List<int> { 4, 3, 3, 7, 7, 7, 7, 3, 3, 12, 8, 4, 4, 13, 7, 7, 7, 12, 12, 3, 3, 3, 2, 2, 2, 2, 2, 2, 11, 11, 11, 11, 11, 11, 4, 4, 4, 4, 4, 4, 10, 10, 13, 13, 5, 5, 5, 5, 15, 15, 15, 15, 3, 3, 3, 15, 15, 15, 3, 3, 3, 16, 6, 2, 2, 18, 3, 2, 2, 2, 2, 2, 9, 9, 9, 3, 2, 2, 2, 2, 2, 3, 3, 12, 6, 3, };

            foreach (PeptideSpectralMatch psm in allPsmsArray)
            {
                if (psm != null)
                {
                    foreach (var (Notch, Peptide) in psm.BestMatchingPeptides)
                    {
                        longestSeriesObserved.Add(PeptideSpectralMatch.GetLongestIonSeriesBidirectional(psm.PeptidesToMatchingFragments, Peptide));
                    }
                }
            }

            Assert.IsTrue(longestSeriesExpected.SequenceEqual(longestSeriesObserved));
        }

        [Test]
        public static void TestDecoyContaminantsFilter()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPSMOutput");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");

            //Filter decoys
            SearchTask searchTaskDecoy = new SearchTask();

            searchTaskDecoy.SearchParameters.WriteDecoys = false;

            var engineDecoy = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("DecoyTest", searchTaskDecoy) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engineDecoy.Run();

            string psmFileDecoy = Path.Combine(outputFolder, @"DecoyTest\AllPSMs.psmtsv");
            var linesDecoy = File.ReadAllLines(psmFileDecoy);
            Assert.That(linesDecoy.Length == 8);

            //Filter contaminants
            SearchTask searchTaskContaminant = new SearchTask();

            searchTaskContaminant.SearchParameters.WriteContaminants = false;

            var engineContaminant = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("ContaminantTest", searchTaskContaminant) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engineContaminant.Run();

            string psmFileContaminant = Path.Combine(outputFolder, @"ContaminantTest\AllPSMs.psmtsv");
            var linesContaminant = File.ReadAllLines(psmFileContaminant);
            Assert.That(linesContaminant.Length == 11);

            string proteinFileContaminant = Path.Combine(outputFolder, @"ContaminantTest\AllQuantifiedProteinGroups.tsv");
            var linesContaminantProtein = File.ReadAllLines(proteinFileContaminant);
            Assert.That(linesContaminantProtein.Length == 7);

            var engineContaminant2 = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("ContaminantTest", searchTaskContaminant) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, true) }, outputFolder);
            engineContaminant2.Run();

            var linesContaminant2 = File.ReadAllLines(psmFileContaminant);
            Assert.That(linesContaminant2.Length == 1);

            var linesContaminantProtein2 = File.ReadAllLines(proteinFileContaminant);
            Assert.That(linesContaminantProtein2.Length == 1);

            //Filter contaminants and decoys
            SearchTask searchTaskDecoyContaminant = new SearchTask();
            SearchTask searchTaskDecoyContaminant2 = new SearchTask();

            searchTaskDecoyContaminant2.SearchParameters.WriteContaminants = false;
            searchTaskDecoyContaminant2.SearchParameters.WriteDecoys = false;

            var engineDecoyContaminant = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("DecoyContaminantTest", searchTaskDecoyContaminant) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, true) }, outputFolder);
            engineContaminant.Run();

            string psmFileDecoyContaminant = Path.Combine(outputFolder, @"DecoyContaminantTest\AllPSMs.psmtsv");
            var linesDecoyContaminant = File.ReadAllLines(psmFileContaminant);
            Assert.That(linesContaminant.Length == 11);

            var engineDecoyContaminant2 = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("DecoyContaminantTest", searchTaskDecoyContaminant2) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, true) }, outputFolder);
            engineContaminant2.Run();

            string psmFileDecoyContaminant2 = Path.Combine(outputFolder, @"DecoyContaminantTest\AllPSMs.psmtsv");
            var linesDecoyContaminant2 = File.ReadAllLines(psmFileContaminant);
            Assert.That(linesContaminant2.Length == 1);

            //No filter
            SearchTask searchTask = new SearchTask();

            var engine = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("NoFilterTest", searchTask) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, true) }, outputFolder);
            engine.Run();

            string psmFile = Path.Combine(outputFolder, @"NoFilterTest\AllPSMs.psmtsv");
            var lines = File.ReadAllLines(psmFile);
            Assert.That(lines.Length == 11);
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestPsmMatchingToTargetAndDecoyWithSameSequence()
        {
            List<Modification> mods = new List<Modification>();
            CommonParameters commonParameters = new CommonParameters();

            PeptideWithSetModifications target = new Protein("PEPTIDE", "TARGET").Digest(commonParameters.DigestionParams, mods, mods).First();
            PeptideWithSetModifications decoy = new Protein("PEPTIDE", "DECOY", isDecoy: true).Digest(commonParameters.DigestionParams, mods, mods).First();

            MsDataFile msDataFile = new TestDataFile(target);
            MsDataScan msDataScan = msDataFile.GetOneBasedScan(2);
            Ms2ScanWithSpecificMass scanWithMass = new Ms2ScanWithSpecificMass(msDataScan, 4, 1, null, commonParameters);

            PeptideSpectralMatch psm = new PeptideSpectralMatch(target, 0, 1, 1, scanWithMass, commonParameters, null);
            psm.AddOrReplace(decoy, 1, 0, true, null, 0);

            Assert.AreEqual(2, psm.BestMatchingPeptides.Count());
            Assert.That(psm.BestMatchingPeptides.Any(p => p.Peptide.Protein.IsDecoy));

            psm.ResolveAllAmbiguities();

            Assert.AreEqual(1, psm.BestMatchingPeptides.Count());
            Assert.That(psm.BestMatchingPeptides.All(p => !p.Peptide.Protein.IsDecoy));
            Assert.That(!psm.IsDecoy);
        }

        [Test]
        public static void TestPsmMatchingToTargetAndDecoyWithDifferentSequences()
        {
            CommonParameters commonParameters = new CommonParameters();
            List<Modification> mods = new List<Modification>();

            PeptideWithSetModifications target = new Protein("PEPTIDE", "").Digest(commonParameters.DigestionParams, mods, mods).First();
            PeptideWithSetModifications decoy = new Protein("PEPTIDEL", "", isDecoy: true).Digest(commonParameters.DigestionParams, mods, mods).First();

            MsDataFile msDataFile = new TestDataFile(target);
            MsDataScan msDataScan = msDataFile.GetOneBasedScan(2);
            Ms2ScanWithSpecificMass scanWithMass = new Ms2ScanWithSpecificMass(msDataScan, 4, 1, null, commonParameters);

            PeptideSpectralMatch psm = new PeptideSpectralMatch(target, 0, 1, 1, scanWithMass, commonParameters, null);
            psm.AddOrReplace(decoy, 1, 0, true, null, 0);

            Assert.AreEqual(2, psm.BestMatchingPeptides.Count());
            Assert.That(psm.BestMatchingPeptides.Any(p => p.Peptide.Protein.IsDecoy));

            psm.ResolveAllAmbiguities();

            Assert.AreEqual(2, psm.BestMatchingPeptides.Count());
            Assert.That(psm.IsDecoy);

            List<(string fileName, CommonParameters fileSpecificParameters)> fsp = new List<(string fileName, CommonParameters fileSpecificParameters)> { ("filename", commonParameters) };

            new FdrAnalysisEngine(new List<PeptideSpectralMatch> { psm }, 1, new CommonParameters(), fsp, new List<string>()).Run();
            Assert.AreEqual(0.5, psm.FdrInfo.CumulativeDecoy);
        }

        [Test]
        public static void TestPsmCount()
        {
            Protein p1 = new Protein("PEPTIDE", null);
            CommonParameters commonParameters = new CommonParameters();
            PeptideWithSetModifications pep1 = p1.Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()).ToList().First();

            Protein p2 = new Protein("PEPTIDE", null);
            PeptideWithSetModifications pep2 = p2.Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()).ToList().First();

            Protein p3 = new Protein("PEPTIDE", null);
            PeptideWithSetModifications pep3 = p3.Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()).ToList().First();

            TestDataFile t = new TestDataFile(new List<PeptideWithSetModifications> { pep1, pep2, pep3 });

            MsDataScan mzLibScan1 = t.GetOneBasedScan(2);
            Ms2ScanWithSpecificMass scan1 = new Ms2ScanWithSpecificMass(mzLibScan1, 0, 1, null, new CommonParameters());
            PeptideSpectralMatch psm1 = new PeptideSpectralMatch(pep1, 0, 0, 0, scan1, commonParameters, new List<MatchedFragmentIon>());

            MsDataScan mzLibScan2 = t.GetOneBasedScan(4);
            Ms2ScanWithSpecificMass scan2 = new Ms2ScanWithSpecificMass(mzLibScan2, 0, 1, null, new CommonParameters());
            PeptideSpectralMatch psm2 = new PeptideSpectralMatch(pep2, 0, 0, 0, scan2, commonParameters, new List<MatchedFragmentIon>());

            MsDataScan mzLibScan3 = t.GetOneBasedScan(6);
            Ms2ScanWithSpecificMass scan3 = new Ms2ScanWithSpecificMass(mzLibScan3, 0, 1, null, new CommonParameters());
            PeptideSpectralMatch psm3 = new PeptideSpectralMatch(pep3, 0, 0, 0, scan3, commonParameters, new List<MatchedFragmentIon>());

            psm1.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0); // valid psm
            psm1.ResolveAllAmbiguities();

            psm2.SetFdrValues(0, 0, 0.02, 0, 0, 0, 0, 0); // psm above fdr cutoff
            psm2.ResolveAllAmbiguities();

            psm3.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0); // ambiguous psm

            var allPsms = new List<PeptideSpectralMatch> { psm1, psm2, psm3 };

            List<(string fileName, CommonParameters fileSpecificParameters)> fsp = new List<(string fileName, CommonParameters fileSpecificParameters)> { ("filename", new CommonParameters()) };

            var fdrEngine = new FdrAnalysisEngine(allPsms, 0, new CommonParameters(), fsp, new List<string>());

            fdrEngine.CountPsm();
            var psmGroups = allPsms.Where(psm => psm.FullSequence != null && psm.PsmCount > 0).GroupBy(p => p.FullSequence).ToList();
            Assert.That(psmGroups.First().Count() == 2);
            Assert.That(psmGroups.First().First().PsmCount == 1);

            psm2.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
            psm3.ResolveAllAmbiguities();

            fdrEngine.CountPsm();
            psmGroups = allPsms.Where(psm => psm.FullSequence != null && psm.PsmCount > 0).GroupBy(p => p.FullSequence).ToList();
            Assert.That(psmGroups.First().Count() == 3);
        }

        [Test]
        public static void TestPsmCount2()
        {
            SearchTask task = new SearchTask();
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPsmCount2");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            Directory.CreateDirectory(outputFolder);

            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { myFile }, "test");

            var peptides = File.ReadAllLines(Path.Combine(outputFolder, @"AllPeptides.psmtsv"));
            var header = peptides[0].Split(new char[] { '\t' }).ToArray();
            int indexOfPsmCountInTsv = Array.IndexOf(header, PsmTsvHeader.PsmCount);
            int indexOfQValueInTsv = Array.IndexOf(header, PsmTsvHeader.QValue);
            Assert.That(indexOfPsmCountInTsv >= 0);
            Assert.That(indexOfQValueInTsv >= 0);

            var psmsFromTsv = PsmTsvReader.ReadTsv(Path.Combine(outputFolder, @"AllPSMs.psmtsv"), out var warnings);
            var psmsGroupedBySequence = psmsFromTsv.GroupBy(p => p.FullSequence).ToList();
            Assert.AreEqual(psmsGroupedBySequence.Count, peptides.Length - 1);

            for (int i = 0; i < psmsGroupedBySequence.Count; i++)
            {
                var peptideLine = peptides[i + 1];
                var split = peptideLine.Split(new char[] { '\t' });

                int psmCount = psmsGroupedBySequence[i].Count(p => p.QValue <= 0.01);
                int psmCountWrittenToPeptidesFile = int.Parse(split[indexOfPsmCountInTsv]);

                Assert.AreEqual(psmCount, psmCountWrittenToPeptidesFile);
            }

            Directory.Delete(outputFolder, true);
            Directory.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"Task Settings"), true);
        }

        [Test]
        public static void PsmtsvTest()
        {
            Type type = typeof(PsmFromTsv);
            PropertyInfo[] properties = type.GetProperties();
        }

        [Test]
        public static void TestPsmAddOrReplace()
        {
            Ms2ScanWithSpecificMass scanB = new Ms2ScanWithSpecificMass(
                new MsDataScan(
                    new MzSpectrum(new double[] { }, new double[] { }, false),
                    2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null),
                100, 1, null, new CommonParameters(), null);

            PeptideSpectralMatch psm1 = new PeptideSpectralMatch(new PeptideWithSetModifications(new Protein("PEPTIDE", "ACCESSION", "ORGANISM"), new DigestionParams(), 1, 2, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0), 0, 10, 1, scanB, new CommonParameters(), new List<MatchedFragmentIon>(), 0);

            PeptideWithSetModifications pwsm = new PeptideWithSetModifications(new Protein("PEPTIDE", "ACCESSION", "ORGANISM"), new DigestionParams(), 1, 2, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

            psm1.AddOrReplace(pwsm, 11, 1, true, new List<MatchedFragmentIon>(), 0);

            Assert.AreEqual(1, psm1.BestMatchingPeptides.Count());

            Assert.AreEqual(11, psm1.Score);

            Assert.AreEqual(10, psm1.RunnerUpScore);
            Assert.AreEqual(1, psm1.DeltaScore);
        }

        [Test]
        public static void TestComplementaryIons()
        {
            Ms2ScanWithSpecificMass scanB = new Ms2ScanWithSpecificMass(
                new MsDataScan(
                    new MzSpectrum(new double[] { }, new double[] { }, false),
                    2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null),
                100, 1, null, new CommonParameters(), null);

            PeptideSpectralMatch psm1 = new PeptideSpectralMatch(new PeptideWithSetModifications(new Protein("PEPTIDE", "ACCESSION", "ORGANISM"), new DigestionParams(), 1, 2, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0), 0, 10, 1, scanB, new CommonParameters(), new List<MatchedFragmentIon>(), 0);

            PeptideWithSetModifications pwsm = new PeptideWithSetModifications(new Protein("PEPTIDE", "ACCESSION", "ORGANISM"), new DigestionParams(), 1, 2, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

            int count = PeptideSpectralMatch.GetCountComplementaryIons(psm1.PeptidesToMatchingFragments, pwsm);

            //No Matched Fragment Ions Returns 0
            Assert.AreEqual(0, count);

            count = PeptideSpectralMatch.GetCountComplementaryIons(null, pwsm);
            //PeptidesToMatchingFragments Null Returns 0
            Assert.AreEqual(0, count);

            List<Product> myProducts = new List<Product>();
            pwsm.Fragment(DissociationType.HCD, FragmentationTerminus.Both, myProducts);
            List<MatchedFragmentIon> mfiList = new List<MatchedFragmentIon>();
            //foreach (Product prod in myProducts)
            for (int i = 0; i < myProducts.Count; i++)
            {
                var prod = myProducts[i];
                mfiList.Add(new MatchedFragmentIon(ref prod, 1, 1, 1));
            }

            Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>> PTMF = new Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>();
            PTMF.Add(pwsm, mfiList);

            count = PeptideSpectralMatch.GetCountComplementaryIons(PTMF, pwsm);
            //PeptidesToMatchingFragments Contains one N and one C ion so intersection Returns 1
            Assert.AreEqual(1, count);
        }

        [Test]
        public static void Test_PSM_GetLongestIonSeries_NullChecks()
        {
            Ms2ScanWithSpecificMass scanB = new Ms2ScanWithSpecificMass(
                new MsDataScan(
                    new MzSpectrum(new double[] { }, new double[] { }, false),
                    2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null),
                100, 1, null, new CommonParameters(), null);

            int longestSeries = 0;

            PeptideWithSetModifications pwsm = new PeptideWithSetModifications(new Protein("PEPTIDE", "ACCESSION", "ORGANISM"), new DigestionParams(), 1, 2, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

            //PeptidesToMatchingFragments == null returns 1
            longestSeries = PeptideSpectralMatch.GetLongestIonSeriesBidirectional(null, pwsm);
            Assert.AreEqual(1, longestSeries);

            //matchedFragments == null returns 1
            Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>> PeptidesToMatchingFragments = new Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>();
            PeptidesToMatchingFragments.Add(pwsm, null);

            longestSeries = PeptideSpectralMatch.GetLongestIonSeriesBidirectional(PeptidesToMatchingFragments, pwsm);
            Assert.AreEqual(1, longestSeries);
        }
    }
}