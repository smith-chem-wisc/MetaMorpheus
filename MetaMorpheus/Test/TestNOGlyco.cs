using EngineLayer;
using EngineLayer.GlycoSearch;
using MassSpectrometry;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public class TestNOGlyco
    {
        [Test]
        public static void NandO_GlycoSearchIndividualFileFolderOutputTest()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            string proteinDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\N_O_glycoWithFileSpecific\\FourMucins_NoSigPeps_FASTA.fasta");
            string spectraFileDirctory = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\N_O_glycoWithFileSpecific");
            List<string> rawFilePaths = Directory.GetFiles(spectraFileDirctory).Where(p => p.Contains("mzML")).ToList();

            // run task
            CommonParameters commonParameters = new(dissociationType: DissociationType.HCD, ms2childScanDissociationType: DissociationType.EThcD);

            Directory.CreateDirectory(outputFolder);
            var glycoSearchTask = new GlycoSearchTask()
            {
                CommonParameters = commonParameters,
                _glycoSearchParameters = new GlycoSearchParameters()
                {
                    OGlycanDatabasefile = "OGlycan.gdb",
                    NGlycanDatabasefile = "NGlycan_prunedForTesting.gdb",
                    GlycoSearchType = GlycoSearchType.N_O_GlycanSearch,
                    OxoniumIonFilt = true,
                    DecoyType = DecoyType.Reverse,
                    GlycoSearchTopNum = 50,
                    MaximumOGlycanAllowed = 3,
                    DoParsimony = true,
                    WriteContaminants = true,
                    WriteDecoys = true,
                    WriteIndividualFiles = true
                }
            };
            GlobalVariables.NGlycanDatabasePaths.Add(Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", @"NGlycan_prunedForTesting.gdb"));
            glycoSearchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(proteinDatabase, false) }, rawFilePaths, "");
            Glycan[] globalO = GlycanBox.GlobalOGlycans;
            var globalN = GlycanBox.GlobalNGlycans;


            Directory.Delete(outputFolder, true);
        }


        [Test]
        public static void TestNOGlycanBox()
        {
            // Loading the glycan databases
            string oglycanPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", @"OGlycan.gdb");
            string nglycanPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", @"NGlycan_prunedForTesting.gdb");

            // reading the O-Glycan database and storing it in GlycanBox.GlobalOGlycans
            GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(oglycanPath, true, true).ToArray();
            Assert.That(GlycanBox.GlobalOGlycans.Count() == 24);

            // reading the N-Glycan database and storing it in GlycanBox.GlobalNGlycans with negative keys
            GlycanBox.GlobalNGlycans = new Dictionary<int, Glycan>();
            var nGlycans = GlycanDatabase.LoadGlycan(nglycanPath, true, true).ToArray();
            int indexForNGlycan = -1;
            foreach (var nGlycan in nGlycans)
            {
                GlycanBox.GlobalNGlycans.Add(indexForNGlycan, nGlycan);
                indexForNGlycan--;
            }
            Assert.That(GlycanBox.GlobalNGlycans.Count == 20);
            Assert.That(!GlycanBox.GlobalNGlycans.Any(p=>p.Key > 0));

            // Building the NO_glycan boxes
            int _maxGlycanNum = 3;
            GlycanBox.NOGlycanBoxes = GlycanBox.BuildNOGlycanBoxes(_maxGlycanNum, false).OrderBy(p => p.Mass).ToArray();
            
            // Check the number of glycan in each box is less than or equal to max glycan number
            Assert.That(!GlycanBox.NOGlycanBoxes.Any(p=>p.NumberOfMods > 3));
            // Check there is at most one N-glycan in each box (p is negative means Nglycan)
            Assert.That(!GlycanBox.NOGlycanBoxes.Any(p => p.ModIds.Count(p => p < 0) > 1));

            // Check the mass of each box is equal to the sum of the mass of glycans in the box
            foreach (var glycanBox in GlycanBox.NOGlycanBoxes)
            {
                double totalOglycanMass = glycanBox.ModIds.Where(p => p >= 0).Select(p => GlycanBox.GlobalOGlycans[p].Mass / 1E5).Sum();
                double totalNglycanMass = glycanBox.ModIds.Where(p => p < 0).Select(p => GlycanBox.GlobalNGlycans[p].Mass / 1E5).Sum();
                Assert.That(glycanBox.Mass, Is.EqualTo(totalNglycanMass + totalOglycanMass).Within(1E-4));
            }
        }
    }
}
