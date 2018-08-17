using EngineLayer;
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
                full_name: "fullName",
                sequenceVariations: new List<SequenceVariation> { new SequenceVariation(2, "P", "Q", "changed this sequence") })
                    .Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).First();
            MsDataFile myMsDataFile = new TestDataFile(pepWithSetMods, "quadratic");
            MsDataScan scann = myMsDataFile.GetOneBasedScan(2);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(scann, 4, 1, null);
            PeptideSpectralMatch psm = new PeptideSpectralMatch(pepWithSetMods.CompactPeptide(TerminusType.None), 1, 2, 3, scan, digestionParams);

            var t = psm.ToString();
            var tabsepheader = PeptideSpectralMatch.GetTabSeparatedHeader();
            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));

            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>
            {
                { pepWithSetMods.CompactPeptide(TerminusType.None), new HashSet<PeptideWithSetModifications>{ pepWithSetMods } }
            };

            psm.MatchToProteinLinkedPeptides(matching);

            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));

            Tolerance fragmentTolerance = new PpmTolerance(10);
            List<ProductType> lp = new List<ProductType> { ProductType.B };

            new LocalizationEngine(new List<PeptideSpectralMatch> { psm }, lp, myMsDataFile, new CommonParameters(productMassTolerance: fragmentTolerance), new List<string>()).Run();

            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));

            psm.SetFdrValues(6, 6, 6, 6, 6, 6, 0, 0, 0, true);

            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));
        }

        [Test]
        public static void TestPSMOutput()
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

            var engine = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("search", searchTask) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engine.Run();
            
            string psmFile = Path.Combine(outputFolder, @"search\AllPSMs.psmtsv");
            var lines = File.ReadAllLines(psmFile);
            Assert.That(lines.Length == 12);

            var engine2 = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("search", searchTask2) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engine2.Run();

            var lines2 = File.ReadAllLines(psmFile);
            Assert.That(lines2.Length == 7);
        }
    }
}