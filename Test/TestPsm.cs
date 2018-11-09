using EngineLayer;
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
            DigestionParams digestionParams = new DigestionParams();
            PeptideWithSetModifications pepWithSetMods = new Protein(
                "MQQQQQQQ",
                "accession1",
                "org",
                new List<Tuple<string, string>> { new Tuple<string, string>("geneNameType", "geneName") },
                new Dictionary<int, List<Modification>> { { 2, new List<Modification> { new Modification("mod", "mod") } } },
                name: "name",
                fullName: "fullName",
                sequenceVariations: new List<SequenceVariation> { new SequenceVariation(2, "P", "Q", "changed this sequence") })
                    .Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();
            MsDataFile myMsDataFile = new TestDataFile(pepWithSetMods, "quadratic");
            MsDataScan scann = myMsDataFile.GetOneBasedScan(2);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(scann, 4, 1, null, new CommonParameters());

            var theoreticalIons = pepWithSetMods.Fragment(DissociationType.HCD, FragmentationTerminus.Both).ToList();
            var matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan, theoreticalIons, new CommonParameters());
            PeptideSpectralMatch psm = new PeptideSpectralMatch(pepWithSetMods, 1, 2, 3, scan, digestionParams, matchedIons);
            psm.ResolveAllAmbiguities();

            var t = psm.ToString();
            var tabsepheader = PeptideSpectralMatch.GetTabSeparatedHeader();
            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));

            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));

            Tolerance fragmentTolerance = new PpmTolerance(10);
            new LocalizationEngine(new List<PeptideSpectralMatch> { psm }, myMsDataFile, new CommonParameters(productMassTolerance: fragmentTolerance), new List<string>()).Run();

            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));

            psm.SetFdrValues(6, 6, 6, 6, 6, 6, 0, 0, 0, true);

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
            Assert.That(lines.Length == 12);

            var engine2 = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("QValueTest", searchTask2) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engine2.Run();

            var lines2 = File.ReadAllLines(psmFile);
            Assert.That(lines2.Length == 7);
            Directory.Delete(outputFolder, true);
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
            Assert.That(linesDecoy.Length == 9);

            //Filter contaminants
            SearchTask searchTaskContaminant = new SearchTask();

            searchTaskContaminant.SearchParameters.WriteContaminants = false;

            var engineContaminant = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("ContaminantTest", searchTaskContaminant) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engineContaminant.Run();

            string psmFileContaminant = Path.Combine(outputFolder, @"ContaminantTest\AllPSMs.psmtsv");
            var linesContaminant = File.ReadAllLines(psmFileContaminant);
            Assert.That(linesContaminant.Length == 12);

            string proteinFileContaminant = Path.Combine(outputFolder, @"ContaminantTest\AllProteinGroups.tsv");
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
            Assert.That(linesContaminant.Length == 12);

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
            Assert.That(lines.Length == 12);
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestPsmMatchingToTargetAndDecoyWithSameSequence()
        {
            DigestionParams digest = new DigestionParams();
            List<Modification> mods = new List<Modification>();

            PeptideWithSetModifications target = new Protein("PEPTIDE", "").Digest(digest, mods, mods).First();
            PeptideWithSetModifications decoy = new Protein("PEPTIDE", "", isDecoy: true).Digest(digest, mods, mods).First();

            MsDataFile msDataFile = new TestDataFile(target);
            MsDataScan msDataScan = msDataFile.GetOneBasedScan(2);
            Ms2ScanWithSpecificMass scanWithMass = new Ms2ScanWithSpecificMass(msDataScan, 4, 1, null, new CommonParameters());

            PeptideSpectralMatch psm = new PeptideSpectralMatch(target, 0, 1, 1, scanWithMass, digest, null);
            psm.AddOrReplace(decoy, 1, 0, true, null);

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
            DigestionParams digest = new DigestionParams();
            List<Modification> mods = new List<Modification>();

            PeptideWithSetModifications target = new Protein("PEPTIDE", "").Digest(digest, mods, mods).First();
            PeptideWithSetModifications decoy = new Protein("PEPTIDEL", "", isDecoy: true).Digest(digest, mods, mods).First();

            MsDataFile msDataFile = new TestDataFile(target);
            MsDataScan msDataScan = msDataFile.GetOneBasedScan(2);
            Ms2ScanWithSpecificMass scanWithMass = new Ms2ScanWithSpecificMass(msDataScan, 4, 1, null, new CommonParameters());

            PeptideSpectralMatch psm = new PeptideSpectralMatch(target, 0, 1, 1, scanWithMass, digest, null);
            psm.AddOrReplace(decoy, 1, 0, true, null);

            Assert.AreEqual(2, psm.BestMatchingPeptides.Count());
            Assert.That(psm.BestMatchingPeptides.Any(p => p.Peptide.Protein.IsDecoy));

            psm.ResolveAllAmbiguities();

            Assert.AreEqual(2, psm.BestMatchingPeptides.Count());
            Assert.That(psm.IsDecoy);

            new FdrAnalysisEngine(new List<PeptideSpectralMatch> { psm }, 1, new CommonParameters(), new List<string>()).Run();
            Assert.AreEqual(0.5, psm.FdrInfo.CumulativeDecoy);
        }
    }
}
