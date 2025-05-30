using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.Localization;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using TaskLayer;
using UsefulProteomicsDatabases;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Easy.Common.Extensions;
using Omics.BioPolymer;
using Readers;
using static Nett.TomlObjectFactory;
using SpectrumMatchFromTsvHeader = EngineLayer.SpectrumMatchFromTsvHeader;

namespace Test
{
    [TestFixture]
    public static class TestPsm
    {

        [Test]
        public static void TestPsmHeader()
        {
            CommonParameters commonParameters = new CommonParameters();
            var pepWithSetMods = new Protein(
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
            SpectralMatch psm = new PeptideSpectralMatch(pepWithSetMods, 1, 2, 3, scan, commonParameters, matchedIons);
            psm.ResolveAllAmbiguities();

            var t = psm.ToString();
            var tabsepheader = SpectralMatch.GetTabSeparatedHeader();
            Assert.That(psm.ToString().Count(f => f == '\t'), Is.EqualTo(SpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t')));

            Assert.That(psm.ToString().Count(f => f == '\t'), Is.EqualTo(SpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t')));

            Tolerance fragmentTolerance = new PpmTolerance(10);
            new LocalizationEngine(new List<SpectralMatch> { psm }, myMsDataFile, new CommonParameters(productMassTolerance: fragmentTolerance), null, new List<string>()).Run();

            Assert.That(psm.ToString().Count(f => f == '\t'), Is.EqualTo(SpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t')));

            psm.SetFdrValues(6, 6, 6, 6, 6, 0, 0, 0);

            Assert.That(psm.ToString().Count(f => f == '\t'), Is.EqualTo(SpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t')));
        }

        [Test]
        public static void TestResolvedAmbiguitiesWithDiagnosticIon()
        {
            CommonParameters commonParameters = new CommonParameters();

            Protein protein = new Protein("", "accession");
            Modification oxP = new Modification(_originalId: "Oxidation on P", _monoisotopicMass: 16);
            Modification oxI = new Modification(_originalId: "Oxidation on I", _monoisotopicMass: 16);
            PeptideWithSetModifications pepWithOxidationOnP= new("PEP[Variable:Oxidation]TIDEK", 
                new Dictionary<string, Modification> { { "Oxidation", oxP} }, p: protein);
            PeptideWithSetModifications pepWithOxidationOnI = new("PEPTI[Variable:Oxidation]DEK",
                new Dictionary<string, Modification> { { "Oxidation", oxI} }, p: protein);

            MsDataFile myMsDataFile = new TestDataFile();
            MsDataScan scann = myMsDataFile.GetOneBasedScan(1);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(scann, 4, 1, null, new CommonParameters());

            var bIon = new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0), 100, 100, 1);
            var dIon = new MatchedFragmentIon(new Product(ProductType.D, FragmentationTerminus.None, 1, 1, 1, 0), 50, 100, 1);
            var twoIonList = new List<MatchedFragmentIon> { bIon, dIon };

            // Start with a psm matched to a peptide with oxidation on I and one matched product ion
            SpectralMatch psm = new PeptideSpectralMatch(pepWithOxidationOnI, 1, 2, 3, scan, commonParameters, new List<MatchedFragmentIon> { bIon });
            // add a pwsm matched to a peptide with oxidation on P, one matched product ion, and one matched diagnostic ion (scores are identical as diagnostic ions aren't scored)
            psm.AddOrReplace(pepWithOxidationOnP, 2, 1, true, twoIonList);

            Assert.That(psm.BestMatchingBioPolymersWithSetMods.Count(), Is.EqualTo(2));

            psm.ResolveAllAmbiguities();

            // Check that the ion series with the diagnostic ion is set to the MatchedFragmentIons property
            Assert.That(psm.MatchedFragmentIons, Is.EqualTo(twoIonList));
            // Check that the pwsm with the diagnostic ion is first in the BestMatchingBioPolymersWithSetMods list
            Assert.That(psm.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer, Is.EqualTo(pepWithOxidationOnP));
        }

        [Test]
        public static void TestQValueFilter()
        {
            // Output filtering is not performed using default settings for CommonParameters
            SearchTask searchTask = new SearchTask();

            SearchTask searchTask2 = new SearchTask()
            {
                CommonParameters = new CommonParameters
                (
                    qValueThreshold: 0
                ),
                SearchParameters = new SearchParameters
                {
                    WriteHighQValuePsms = false
                }
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
        public static void TestPpmAndDaMassErrors()
        {
            
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_HeLa_04_subset_longestSeq.mzML");
            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters CommonParameters = new CommonParameters();
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
            fsp.Add((origDataFile, CommonParameters));
            var myMsDataFile = myFileManager.LoadFile(origDataFile, CommonParameters);
            var searchModes = new SinglePpmAroundZeroSearchMode(5);
            List<Protein> proteinList = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\hela_snip_for_unitTest.fasta"), true, DecoyType.Reverse, false, out var dbErrors,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex, -1);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();
            SpectralMatch[] allPsmsArray = new SpectralMatch[listOfSortedms2Scans.Length];
            bool writeSpetralLibrary = false;
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null, proteinList, searchModes,
                new CommonParameters(), null, null, new List<string>(), writeSpetralLibrary).Run();
            FdrAnalysisResults fdrResultsClassicDelta = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArray.Where(p => p != null).ToList(), 1,
                CommonParameters, fsp, new List<string>()).Run());
            var nonNullPsms = allPsmsArray.Where(p => p != null).ToList();

            foreach (PeptideSpectralMatch psm in nonNullPsms)
            {
                double daError =
                    Math.Round(psm.ScanPrecursorMass - psm.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer.MonoisotopicMass, 5);
                Assert.That(psm.PrecursorMassErrorDa.First(), Is.EqualTo(daError).Within(0.01));

                double ppmError = Math.Round((psm.ScanPrecursorMass - psm.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer.MonoisotopicMass) / psm.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer.MonoisotopicMass * 1e6, 5);
                Assert.That(psm.PrecursorMassErrorPpm.First(), Is.EqualTo(ppmError).Within(0.1));
            }
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
            SpectralMatch[] allPsmsArray = new SpectralMatch[listOfSortedms2Scans.Length];
            bool writeSpetralLibrary = false;
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null, proteinList, searchModes,
                new CommonParameters(), null, null, new List<string>(), writeSpetralLibrary).Run();

            List<int> longestSeriesObserved = new List<int>();
            List<int> longestSeriesExpected = new List<int> { 4, 3, 3, 7, 7, 7, 7, 3, 3, 12, 8, 4, 4, 13, 7, 7, 7, 12, 12, 3, 3, 3, 2, 2, 2, 2, 2, 2, 11, 11, 11, 11, 11, 11, 4, 4, 4, 4, 4, 4, 10, 10, 13, 13, 5, 5, 5, 5, 15, 15, 15, 15, 3, 3, 3, 15, 15, 15, 3, 3, 3, 16, 6, 2, 2, 18, 3, 2, 2, 2, 2, 2, 9, 9, 9, 3, 2, 2, 2, 2, 2, 3, 3, 12, 6, 3, };

            foreach (SpectralMatch psm in allPsmsArray)
            {
                if (psm != null)
                {
                    foreach (var bestMatch in psm.BestMatchingBioPolymersWithSetMods)
                    {
                        longestSeriesObserved.Add(SpectralMatch.GetLongestIonSeriesBidirectional(bestMatch));
                    }
                }
            }
            Assert.That(longestSeriesExpected.SequenceEqual(longestSeriesObserved));
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

            var target = new Protein("PEPTIDE", "TARGET").Digest(commonParameters.DigestionParams, mods, mods).First();
            var decoy = new Protein("PEPTIDE", "DECOY", isDecoy: true).Digest(commonParameters.DigestionParams, mods, mods).First();

            MsDataFile msDataFile = new TestDataFile(target);
            MsDataScan msDataScan = msDataFile.GetOneBasedScan(2);
            Ms2ScanWithSpecificMass scanWithMass = new Ms2ScanWithSpecificMass(msDataScan, 4, 1, null, commonParameters);

            SpectralMatch psm = new PeptideSpectralMatch(target, 0, 1, 1, scanWithMass, commonParameters, null);
            psm.AddOrReplace(decoy, 1, 0, true, null);

            Assert.That(psm.BestMatchingBioPolymersWithSetMods.Count(), Is.EqualTo(2));
            Assert.That(psm.BestMatchingBioPolymersWithSetMods.Any(p => p.SpecificBioPolymer.Parent.IsDecoy));

            psm.ResolveAllAmbiguities();

            Assert.That(psm.BestMatchingBioPolymersWithSetMods.Count(), Is.EqualTo(1));
            Assert.That(psm.BestMatchingBioPolymersWithSetMods.All(p => !p.SpecificBioPolymer.Parent.IsDecoy));
            Assert.That(!psm.IsDecoy);
        }

        [Test]
        public static void TestPsmMatchingToTargetAndDecoyWithDifferentSequences()
        {
            CommonParameters commonParameters = new CommonParameters();
            List<Modification> mods = new List<Modification>();

            var target = new Protein("PEPTIDE", "").Digest(commonParameters.DigestionParams, mods, mods).First();
            var decoy = new Protein("PEPTIDEL", "", isDecoy: true).Digest(commonParameters.DigestionParams, mods, mods).First();

            MsDataFile msDataFile = new TestDataFile(target);
            MsDataScan msDataScan = msDataFile.GetOneBasedScan(2);
            Ms2ScanWithSpecificMass scanWithMass = new Ms2ScanWithSpecificMass(msDataScan, 4, 1, null, commonParameters);

            SpectralMatch psm = new PeptideSpectralMatch(target, 0, 1, 1, scanWithMass, commonParameters, null);
            psm.AddOrReplace(decoy, 1, 0, true, null);

            Assert.That(psm.BestMatchingBioPolymersWithSetMods.Count(), Is.EqualTo(2));
            Assert.That(psm.BestMatchingBioPolymersWithSetMods.Any(p => p.SpecificBioPolymer.Parent.IsDecoy));

            psm.ResolveAllAmbiguities();

            Assert.That(psm.BestMatchingBioPolymersWithSetMods.Count(), Is.EqualTo(2));
            Assert.That(psm.IsDecoy);

            List<(string fileName, CommonParameters fileSpecificParameters)> fsp = new List<(string fileName, CommonParameters fileSpecificParameters)> { ("filename", commonParameters) };

            new FdrAnalysisEngine(new List<SpectralMatch> { psm }, 1, new CommonParameters(), fsp, new List<string>()).Run();
            Assert.That(psm.FdrInfo.CumulativeDecoy, Is.EqualTo(0.5));

        }

        [Test]
        public static void TestPsmCount()
        {
            Protein p1 = new Protein("PEPTIDE", null);
            CommonParameters commonParameters = new CommonParameters();
            var pep1 = p1.Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()).ToList().First();

            Protein p2 = new Protein("PEPTIDE", null);
            var pep2 = p2.Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()).ToList().First();

            Protein p3 = new Protein("PEPTIDE", null);
            var pep3 = p3.Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()).ToList().First();

            TestDataFile t = new TestDataFile(new List<IBioPolymerWithSetMods> { pep1, pep2, pep3 });

            MsDataScan mzLibScan1 = t.GetOneBasedScan(2);
            Ms2ScanWithSpecificMass scan1 = new Ms2ScanWithSpecificMass(mzLibScan1, 0, 1, null, new CommonParameters());
            SpectralMatch psm1 = new PeptideSpectralMatch(pep1, 0, 0, 0, scan1, commonParameters, new List<MatchedFragmentIon>());

            MsDataScan mzLibScan2 = t.GetOneBasedScan(4);
            Ms2ScanWithSpecificMass scan2 = new Ms2ScanWithSpecificMass(mzLibScan2, 0, 1, null, new CommonParameters());
            SpectralMatch psm2 = new PeptideSpectralMatch(pep2, 0, 0, 0, scan2, commonParameters, new List<MatchedFragmentIon>());

            MsDataScan mzLibScan3 = t.GetOneBasedScan(6);
            Ms2ScanWithSpecificMass scan3 = new Ms2ScanWithSpecificMass(mzLibScan3, 0, 1, null, new CommonParameters());
            SpectralMatch psm3 = new PeptideSpectralMatch(pep3, 0, 0, 0, scan3, commonParameters, new List<MatchedFragmentIon>());

            psm1.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0); // valid psm
            psm1.ResolveAllAmbiguities();

            psm2.SetFdrValues(0, 0, 0.02, 0, 0, 0, 0, 0); // psm above fdr cutoff
            psm2.ResolveAllAmbiguities();

            psm3.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0); // ambiguous psm

            var allPsms = new List<SpectralMatch> { psm1, psm2, psm3 };

            List<(string fileName, CommonParameters fileSpecificParameters)> fsp = new List<(string fileName, CommonParameters fileSpecificParameters)> { ("filename", new CommonParameters()) };

            var fdrEngine = new FdrAnalysisEngine(allPsms, 0, new CommonParameters(), fsp, new List<string>());

            fdrEngine.CountPsm(allPsms);
            var psmGroups = allPsms.Where(psm => psm.FullSequence != null && psm.PsmCount > 0).GroupBy(p => p.FullSequence).ToList();
            Assert.That(psmGroups.First().Count() == 2);
            Assert.That(psmGroups.First().First().PsmCount == 1);

            psm2.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
            psm3.ResolveAllAmbiguities();

            fdrEngine.CountPsm(allPsms);
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

            List<string> peptides = File.ReadAllLines(Path.Combine(outputFolder, @"AllPeptides.psmtsv")).ToList();
            var header = peptides[0].Split(new char[] { '\t' }).ToArray();
            int indexOfPsmCountInTsv = Array.IndexOf(header, SpectrumMatchFromTsvHeader.SpectrumMatchCount);
            int indexOfQValueInTsv = Array.IndexOf(header, SpectrumMatchFromTsvHeader.QValue);
            Assert.That(indexOfPsmCountInTsv >= 0);
            Assert.That(indexOfQValueInTsv >= 0);

            peptides.RemoveAt(0);//delete header line
            var allPeptidesQvalueBelowCutoff = 0;
            foreach (var peptide in peptides)
            {
                string qValueString = peptide.Split('\t')[indexOfQValueInTsv].ToString();
                double qValue = Convert.ToDouble(qValueString);
                if (qValue < 0.01)
                {
                    allPeptidesQvalueBelowCutoff++;
                }
            }

            var psmsFromTsv = SpectrumMatchTsvReader.ReadPsmTsv(Path.Combine(outputFolder, @"AllPSMs.psmtsv"), out var warnings);
            var allUnambiguousPsms = psmsFromTsv.Where(psm => psm.FullSequence != null);
            var unambiguousPsmsLessThanOnePercentFdr = allUnambiguousPsms.Where(psm =>
                    psm.QValue <= 0.01)
                .GroupBy(p => p.FullSequence).ToList();
            Assert.That(unambiguousPsmsLessThanOnePercentFdr.Count, Is.EqualTo(allPeptidesQvalueBelowCutoff));

            // Test for precursorIntensity in PsmFromTsv
            int scanNumber = 117;
            var psmOfInterest = psmsFromTsv.First(psm => psm.PrecursorScanNum == scanNumber);
            Assert.That(psmOfInterest.PrecursorIntensity == 4634473.5);

            var lines = File.ReadAllLines(Path.Combine(outputFolder, @"AllPSMs.psmtsv")).ToList();
            int indexOfPrecursorIntensity = Array.IndexOf(header, SpectrumMatchFromTsvHeader.PrecursorIntensity);
            var copy = lines[lines.Count - 1].Split('\t').ToArray();
            copy[indexOfPrecursorIntensity] = copy[indexOfPrecursorIntensity] + "a";
            string line = string.Join("\t", copy);
            lines.Add(line);
            File.WriteAllLines(Path.Combine(outputFolder, @"TestInvalidPSMs.psmtsv"), lines.ToArray());
            var psmsFromTsvInvalid = SpectrumMatchTsvReader.ReadPsmTsv(Path.Combine(outputFolder, @"TestInvalidPSMs.psmtsv"), out var warnings1);
            var psmInvalid = psmsFromTsvInvalid[psmsFromTsvInvalid.Count - 1];
            Assert.That(psmInvalid.PrecursorIntensity, Is.EqualTo(null));

            //Test for precursorIntensity and precursorEnvelopePeakCount in SpectralMatch
            List<Protein> proteinList = ProteinDbLoader.LoadProteinFasta(myDatabase, true, DecoyType.Reverse, false, out List<string> errors);
            var fsp = new List<(string, CommonParameters)>();
            CommonParameters commonParameters = new CommonParameters();
            fsp.Add(("SmallCalibratible_Yeast.mzML", commonParameters));

            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters CommonParameters = new CommonParameters();
            var myMsDataFile = myFileManager.LoadFile(myFile, CommonParameters);
            var arrayOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, myFile, commonParameters).OrderBy(b => b.PrecursorMass).ToArray();
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var searchModes = new SinglePpmAroundZeroSearchMode(5);
            bool writeSpectralLibrary = false;
            SpectralMatch[] allPsmsArray = new SpectralMatch[arrayOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, arrayOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchModes, CommonParameters, null, null, new List<string>(), writeSpectralLibrary).Run();

            List<string> results = File.ReadAllLines(Path.Combine(outputFolder, @"results.txt")).ToList();

            string peptideCountFromResultsString = results.FirstOrDefault(r => r.Contains("All target peptides with q-value <= 0.01: "));
            double peptideCountFromResults = Convert.ToDouble(peptideCountFromResultsString?.Split(':')[1].ToString());
            Assert.That(allPeptidesQvalueBelowCutoff, Is.EqualTo(peptideCountFromResults));
            Directory.Delete(outputFolder, true);
            Directory.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"Task Settings"), true);
        }

        [Test]
        public static void TestPsmAddOrReplace()
        {
            Ms2ScanWithSpecificMass scanB = new Ms2ScanWithSpecificMass(
                new MsDataScan(
                    new MzSpectrum(new double[] { }, new double[] { }, false),
                    2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null),
                100, 1, null, new CommonParameters(), null);

            SpectralMatch psm1 = new PeptideSpectralMatch(new PeptideWithSetModifications(new Protein("PEPTIDE", "ACCESSION", "ORGANISM"), new DigestionParams(), 1, 2, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0), 0, 10, 1, scanB, new CommonParameters(), new List<MatchedFragmentIon>());

            PeptideWithSetModifications pwsm = new PeptideWithSetModifications(new Protein("PEPTIDE", "ACCESSION", "ORGANISM"), new DigestionParams(), 1, 2, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

            psm1.AddOrReplace(pwsm, 11, 1, true, new List<MatchedFragmentIon>());

            Assert.That(psm1.BestMatchingBioPolymersWithSetMods.Count(), Is.EqualTo(1));

            Assert.That(psm1.Score, Is.EqualTo(11));

            Assert.That(psm1.RunnerUpScore, Is.EqualTo(10));
            Assert.That(psm1.DeltaScore, Is.EqualTo(1));
        }

        [Test]
        public static void TestComplementaryIons()
        {
            Ms2ScanWithSpecificMass scanB = new Ms2ScanWithSpecificMass(
                new MsDataScan(
                    new MzSpectrum(new double[] { }, new double[] { }, false),
                    2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null),
                100, 1, null, new CommonParameters(), null);

            SpectralMatch psm1 = new PeptideSpectralMatch(new PeptideWithSetModifications(new Protein("PEPTIDE", "ACCESSION", "ORGANISM"), new DigestionParams(), 1, 2, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0), 0, 10, 1, scanB, new CommonParameters(), new List<MatchedFragmentIon>());

            PeptideWithSetModifications pwsm = new PeptideWithSetModifications(new Protein("PEPTIDE", "ACCESSION", "ORGANISM"), new DigestionParams(), 1, 2, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

            int count = SpectralMatch.GetCountComplementaryIons([], pwsm);

            //No Matched Fragment Ions Returns 0
            Assert.That(count, Is.EqualTo(0));

            count = SpectralMatch.GetCountComplementaryIons(null, pwsm);
            //BioPolymersWithSetModsToMatchingFragments Null Returns 0
            Assert.That(count, Is.EqualTo(0));

            List<Product> myProducts = new List<Product>();
            pwsm.Fragment(DissociationType.HCD, FragmentationTerminus.Both, myProducts);
            List<MatchedFragmentIon> mfiList = new List<MatchedFragmentIon>();
            //foreach (Product prod in myProducts)
            for (int i = 0; i < myProducts.Count; i++)
            {
                var prod = myProducts[i];
                mfiList.Add(new MatchedFragmentIon(prod, 1, 1, 1));
            }

            count = SpectralMatch.GetCountComplementaryIons(mfiList, pwsm);
            //BioPolymersWithSetModsToMatchingFragments Contains one N and one C ion so intersection Returns 1
            Assert.That(count, Is.EqualTo(1));
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

            //BioPolymersWithSetModsToMatchingFragments == null returns 1
            longestSeries = SpectralMatch.GetLongestIonSeriesBidirectional(null, pwsm);
            Assert.That(longestSeries, Is.EqualTo(1));

            //matchedFragments == null returns 1

            longestSeries = SpectralMatch.GetLongestIonSeriesBidirectional(null, pwsm);
            Assert.That(longestSeries, Is.EqualTo(1));
        }

        [Test]
        public static void TestPSMFragmentCoverage()
        {
            CommonParameters commonParameters = new CommonParameters();

            Protein p1 = new Protein("PEPTIDEPEPTIDE", null);
            var pep1 = p1.Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()).ToList().First();

            Protein p2 = new Protein("GGGGGGGGGGGGGGKPEPTIDEPEPTIDE", null);
            var pep2 = p2.Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()).ToList()[1];

            TestDataFile t = new TestDataFile(new List<IBioPolymerWithSetMods> { pep1, pep2 });

            //psm 1 - test first and last amino acid positions, along with one internal Amino Acid position
            Product productC3 = new Product(ProductType.y, FragmentationTerminus.C, 0, 3, 12, 0);
            Product productC4 = new Product(ProductType.y, FragmentationTerminus.C, 0, 4, 11, 0);
            Product productC7 = new Product(ProductType.y, FragmentationTerminus.C, 0, 7, 8, 0);
            Product productC13 = new Product(ProductType.y, FragmentationTerminus.C, 0, 13, 2, 0);
            Product productN3 = new Product(ProductType.b, FragmentationTerminus.N, 0, 3, 3, 0);
            Product productN4 = new Product(ProductType.b, FragmentationTerminus.N, 0, 4, 4, 0);
            Product productN6 = new Product(ProductType.b, FragmentationTerminus.N, 0, 6, 6, 0);
            Product productN8 = new Product(ProductType.b, FragmentationTerminus.N, 0, 8, 8, 0);
            Product productN13 = new Product(ProductType.b, FragmentationTerminus.N, 0, 13, 13, 0);
            MatchedFragmentIon mfiC3 = new MatchedFragmentIon(productC3, 0, 0, 1);
            MatchedFragmentIon mfiC4 = new MatchedFragmentIon(productC4, 0, 0, 1);
            MatchedFragmentIon mfiC7 = new MatchedFragmentIon(productC7, 0, 0, 1);
            MatchedFragmentIon mfiC13 = new MatchedFragmentIon(productC13, 0, 0, 1);
            MatchedFragmentIon mfiN3 = new MatchedFragmentIon(productN3, 0, 0, 1);
            MatchedFragmentIon mfiN4 = new MatchedFragmentIon(productN4, 0, 0, 1);
            MatchedFragmentIon mfiN6 = new MatchedFragmentIon(productN6, 0, 0, 1);
            MatchedFragmentIon mfiN8 = new MatchedFragmentIon(productN8, 0, 0, 1);
            MatchedFragmentIon mfiN13 = new MatchedFragmentIon(productN13, 0, 0, 1);
            List<MatchedFragmentIon> mfis1 = new List<MatchedFragmentIon> { mfiC3, mfiC4, mfiC7, mfiC13, mfiN3, mfiN4, mfiN6, mfiN8, mfiN13 };
            MsDataScan mzLibScan1 = t.GetOneBasedScan(2);
            Ms2ScanWithSpecificMass scan1 = new Ms2ScanWithSpecificMass(mzLibScan1, 0, 1, null, new CommonParameters());
            SpectralMatch psm1 = new PeptideSpectralMatch(pep1, 0, 0, 0, scan1, commonParameters, mfis1);
            psm1.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0); // valid psm
            psm1.ResolveAllAmbiguities();
            psm1.GetAminoAcidCoverage();
            //First amino acid
            Assert.That(psm1.FragmentCoveragePositionInPeptide.Contains(1));
            //sequential N term Frags
            Assert.That(psm1.FragmentCoveragePositionInPeptide.Contains(4));
            //Last amino acid
            Assert.That(psm1.FragmentCoveragePositionInPeptide.Contains(14));
            //Covered from both directions inclusive
            Assert.That(psm1.FragmentCoveragePositionInPeptide.Contains(8));
            //Covered from both directions exclusive
            Assert.That(psm1.FragmentCoveragePositionInPeptide.Contains(7));
            //Sequential C term Frags
            Assert.That(psm1.FragmentCoveragePositionInPeptide.Contains(11));
            //Not coveredRT
            Assert.That(!psm1.FragmentCoveragePositionInPeptide.Contains(5));


            SpectralMatch psm2 = new PeptideSpectralMatch(pep2, 0, 0, 0, scan1, commonParameters, mfis1);
            psm2.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0); // valid psm
            psm2.ResolveAllAmbiguities();
            psm2.GetAminoAcidCoverage();

            //check that fragment coverage positions are the same
            Assert.That(psm1.FragmentCoveragePositionInPeptide.SequenceEqual(psm2.FragmentCoveragePositionInPeptide));
        }

        [Test]
        public static void TestPrecursorIntensity()
        {
            //Test for Ms2WithSpecificMass
            //1: do the deconvolution and use isotopic envelope to find precursor info
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters CommonParameters = new CommonParameters();
            var myMsDataFile = myFileManager.LoadFile(filePath, CommonParameters);

            var scansWithPrecursors = MetaMorpheusTask._GetMs2Scans(myMsDataFile, filePath, CommonParameters);
            var Ms2Scan1 = scansWithPrecursors[17][1];
            Assert.That(Math.Abs(2889051 - Ms2Scan1.PrecursorIntensity) <= 10);
            Assert.That(Ms2Scan1.PrecursorEnvelopePeakCount, Is.EqualTo(2)); //might not be the correct number of peaks but use it for now

            CommonParameters CommonParameters1 = new CommonParameters(useMostAbundantPrecursorIntensity: false);
            var scansWithPrecursors1 = MetaMorpheusTask._GetMs2Scans(myMsDataFile, filePath, CommonParameters1);
            var Ms2Scan1_2 = scansWithPrecursors1[17][1];
            Assert.That(Math.Abs(3405218 - Ms2Scan1_2.PrecursorIntensity) <= 10);

            //just to look at the envelopes, not relavent to the test
            var msNScans = myMsDataFile.GetAllScansList().ToArray();
            var ms2Scan23 = msNScans.Where(p => p.OneBasedScanNumber == 23).First();
            var precursorSpectrum22 = msNScans.Where(p => p.OneBasedScanNumber == 22).First();
            var envelopes = ms2Scan23.GetIsolatedMassesAndCharges(precursorSpectrum22.MassSpectrum, CommonParameters.PrecursorDeconvolutionParameters);

            //2: use scan header (selectedIonMonoisotopicGuessIntensity) to find precursor info
            CommonParameters CommonParameters2 = new CommonParameters(doPrecursorDeconvolution: false, useProvidedPrecursorInfo: true);
            var scansWithPrecursors2 = MetaMorpheusTask._GetMs2Scans(myMsDataFile, filePath, CommonParameters2);
            var Ms2Scan2 = scansWithPrecursors2[17][0];
            Assert.That(Math.Abs(1.14554e7 - Ms2Scan2.PrecursorIntensity) <= 1000);
            Assert.That(Ms2Scan2.PrecursorEnvelopePeakCount, Is.EqualTo(1));

            //3: use scan header (selectedIonIntensity) to find precursor info 
            MzSpectrum spectrum1 = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 0, 1, 2 }, false);
            MzSpectrum spectrum2 = new MzSpectrum(new double[] { 2, 3, 4 }, new double[] { 1000, 2, 4 }, false);
            MsDataScan[] scans = new MsDataScan[2];
            scans[0] = new MsDataScan(spectrum1, 1, 1, true, Polarity.Positive, 1.0, new MzRange(300, 2000), "scan filter", MZAnalyzerType.Unknown, spectrum1.SumOfAllY, null, null, null);
            scans[1] = new MsDataScan(spectrum2, 2, 2, true, Polarity.Positive, 1, new MzRange(300, 2000), "scan filter", MZAnalyzerType.Unknown, spectrum2.SumOfAllY, 1, new double[,] { }, 
                "nativeId", selectedIonMz: 2, selectedIonChargeStateGuess: 1, selectedIonIntensity: 1000, 1, 1, DissociationType.Unknown, null, null, null);

            var testMsDataFile = new GenericMsDataFile(scans, new SourceFile("no nativeID format", "mzML format",
                    null, null, null));
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(testMsDataFile, "mzMLWithZeros.mzML", false);

            var scansWithPrecursors3 = MetaMorpheusTask._GetMs2Scans(testMsDataFile, "mzMLWithZeros.mzML", CommonParameters2);
            var Ms2Scan3 = scansWithPrecursors3[0][0];
            Assert.That(Ms2Scan3.PrecursorIntensity, Is.EqualTo(1000));
            Assert.That(Ms2Scan3.PrecursorEnvelopePeakCount, Is.EqualTo(1));

            //Test for SpectralMatch
            SearchTask task = new SearchTask();
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            List<Protein> proteinList = ProteinDbLoader.LoadProteinFasta(myDatabase, true, DecoyType.Reverse, false, out List<string> errors);
            var fsp = new List<(string, CommonParameters)>();
            CommonParameters commonParameters = new CommonParameters();
            fsp.Add(("SmallCalibratible_Yeast.mzML", commonParameters));
            var arrayOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, filePath, commonParameters).OrderBy(b => b.PrecursorMass).ToArray();
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            variableModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => CommonParameters.ListOfModsVariable.Contains((b.ModificationType, b.IdWithMotif))).ToList();
            fixedModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => CommonParameters.ListOfModsFixed.Contains((b.ModificationType, b.IdWithMotif))).ToList();
            var searchModes = new SinglePpmAroundZeroSearchMode(5);
            bool writeSpectralLibrary = false;
            MassDiffAcceptor massDiffAcceptor = new DotMassDiffAcceptor("1mm", new List<double> { 0, 1.0029 }, new PpmTolerance(5));
            SpectralMatch[] allPsmsArray = new SpectralMatch[arrayOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, arrayOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, massDiffAcceptor, CommonParameters, fsp, null, new List<string>(), writeSpectralLibrary).Run();

            List<SpectralMatch> psms = new List<SpectralMatch>();
            foreach(SpectralMatch psm in allPsmsArray)
            {
                if (psm != null)
                {
                    psms.Add(psm);
                }
            }
            SpectralMatch psmScan23 = psms.ToArray()[33];
            Assert.That(psmScan23.PrecursorScanEnvelopePeakCount, Is.EqualTo(4));
        }
        }
}