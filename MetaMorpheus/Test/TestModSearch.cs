using Easy.Common.Extensions;
using EngineLayer;
using EngineLayer.GlycoSearch;
using FlashLFQ;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection.Metadata;
using System.Runtime.InteropServices.Marshalling;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Animation;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public class TestModSearch
    {
        [Test]
        public static void TestFragment()
        {
            List<ProductType> customIons = new List<ProductType>() { ProductType.c , ProductType.zDot};
            Protein protein = new Protein("DANNTQFQFTSR", "accession");
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 7);
            var pep = protein.Digest(digestionParams, new List<Modification>(), new List<Modification>());
            HashSet<string> motifs = new HashSet<string> { "Nxs", "Nxt" };
            var sites = GlycoSpectralMatch.GetPossibleModSites(pep.Last(), motifs).Select(p => p.Key).ToList();
            byte[] kind = GlycanDatabase.String2Kind("HexNAc(2)Hex(5)");
            Glycan glycan = new Glycan(kind, "Nxs", GlycanType.N_glycan);
            var glycopep = GlycoPeptides.GenerateGlycopeptide(sites[0], pep.Last(), glycan);
            List<Product> products_1 = GlycoPeptides.GetTheoreticalFragments(DissociationType.EThcD, customIons, pep.Last(), glycopep);
            List<Product> products_2 = new List<Product>(); 
            glycopep.Fragment(DissociationType.EThcD, FragmentationTerminus.Both, products_2);
            products_2 = products_2.Where(p => p.ProductType != ProductType.M).ToList();
            Assert.That(products_1.Count == products_2.Count);
            CollectionAssert.AreNotEqual(products_2, products_1);
            products_1 = products_1.Where(p => p.ProductType != ProductType.y && p.ProductType != ProductType.b)
                .ToList();
            products_2 = products_2.Where(p => p.ProductType != ProductType.y && p.ProductType != ProductType.b)
                .ToList();
            bool isSame = true;
            foreach (var product in products_1)
            {
                if (!products_2.Contains(product))
                {
                    isSame = false;
                }
            }
        }

        [Test]
        public static void TestModSearchEngine()
        {
            string outputFolder = "";

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSearchTaskconfigOGlycoTest_Run.toml"), MetaMorpheusTask.tomlConfig);
            glycoSearchTask._glycoSearchParameters.MaximumOGlycanAllowed = 3;
            glycoSearchTask._glycoSearchParameters.GlycoSearchType = GlycoSearchType.ModSearch;
            DbForTask db = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P16150.fasta"), false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) }, new List<string> { spectraFile }, new List<DbForTask> { db }, outputFolder).Run();
        }

        [Test]
        public static void TestModPair_O()
        {
            // set up output directory, MS data, and database
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);
            DbForTask db = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P16150.fasta"), false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf");

            // Do the O-pair search first
            var task_OPair = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSearchTaskconfigOGlycoTest_Run.toml"), MetaMorpheusTask.tomlConfig);
            task_OPair._glycoSearchParameters.MaximumOGlycanAllowed = 3;
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task_OPair) }, new List<string> { spectraFile }, new List<DbForTask> 
                { db }, outputFolder).Run();
            var db_OSearch = ModBox.GlobalModifications;

            // Store O-pair results
            var peptideHearList = File
                .ReadLines(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData", "Task", "AllPSMs.psmtsv"))
                .First().Split('\t');
            List<List<string>> peptideResults_OSearch = new List<List<string>>();
            foreach (var line in File.ReadAllLines(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData", "Task", "AllPSMs.psmtsv")).Skip(1))
            {
                peptideResults_OSearch.Add(line.Split('\t').ToList());
            }
            Directory.Delete(outputFolder, true);


            // Now run the ModPair search
            var task_ModPair = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSearchTaskconfigOGlycoTest_Run.toml"), MetaMorpheusTask.tomlConfig);
            task_ModPair._glycoSearchParameters.GlycoSearchType = GlycoSearchType.ModSearch;
            task_ModPair._glycoSearchParameters.MaximumOGlycanAllowed = 3;
            task_ModPair._glycoSearchParameters.ListOfInterestedMods = ModBox.GlobalModifications // Use the same glycan database as O-pair search
                .Where(b => b.ModificationType.Equals("O-linked glycosylation"))
                .Select(b => (b.ModificationType, b.IdWithMotif))
                .ToList();
            Directory.CreateDirectory(outputFolder);
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task_ModPair) }, new List<string> { spectraFile }, new List<DbForTask>
            {
                db
            }, outputFolder).Run();
            var db_ModSearch = ModBox.GlobalModifications;
            // Read ModSearch results
            List<List<string>> peptideResults_ModSearch = new List<List<string>>();
            foreach (var line in File.ReadAllLines(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData", "Task", "AllPSMs.psmtsv")).Skip(1))
            {
                peptideResults_ModSearch.Add(line.Split('\t').ToList());
            }


            // Compare O-pair and ModPair results
            Assert.That(db_OSearch.Length == db_ModSearch.Length, "O-pair and ModPair databases are not equal"); // Make sure the same database is used for both searches
            Assert.That(peptideResults_OSearch.Count == peptideResults_ModSearch.Count);
            for (int i = 0; i < peptideResults_OSearch.Count; i++) // Compare the results of O-pair and ModPair search
            {
                for (int j = 0; j < peptideResults_OSearch[i].Count; j++)
                {
                    if (j == 32) // The number of the glycoSites
                    {
                        Assert.That(peptideResults_OSearch[i][j] == "8"); // For O pair search, the number of glycoSites is the number of S and T sites
                        Assert.That(peptideResults_ModSearch[i][j] == "6"); // For ModPair search, the number of glycoSites is the number of the mod's motif in the Box. In this case, the motif is S, then there are 6 S sites.
                    }
                    else
                    {
                        Assert.That(peptideResults_OSearch[i][j], Is.EqualTo(peptideResults_ModSearch[i][j]), $"{peptideHearList[j]} compare: {peptideResults_OSearch[i][j]} is not equal to {peptideResults_ModSearch[i][j]}");
                    }
                }
            }
        }


        [Test]
        public static void TestModPair_N()
        {
            // set up output directory, MS data, and database
            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\Q9C0Y4.fasta"), false);
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\yeast_glycan_25170.mgf");
            Directory.CreateDirectory(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"));

            // first, we run the NPair search (O-pair for Nglyco)
            var task_NPair = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\NGlycanSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            task_NPair._glycoSearchParameters.MaximumOGlycanAllowed = 1;
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task_NPair) }, new List<string> { raw }, new List<DbForTask>
            { db }, Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData")).Run();
            var db_NSearch = GlycanDatabase.LoadGlycan(GlobalVariables.NGlycanDatabasePaths.Where(p => System.IO.Path.GetFileName(p) == task_NPair._glycoSearchParameters.NGlycanDatabasefile).First(), true, false).ToList();

            // Store N-pair results
            var peptideHearList = File.ReadLines(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData", "Task", "AllPSMs.psmtsv")).First().Split('\t');
            var peptideResult_NSearch = File.ReadLines(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData", "Task", "AllPSMs.psmtsv")).Skip(1).First().Split('\t'); // Only one psm in this test
            Directory.Delete(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"), true);


            // Now run the ModPair search
            Directory.CreateDirectory(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"));
            var task_ModPair = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\NGlycanSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            task_ModPair._glycoSearchParameters.GlycoSearchType = GlycoSearchType.ModSearch;
            task_ModPair._glycoSearchParameters.MaximumOGlycanAllowed = 1;
            task_ModPair._glycoSearchParameters.ListOfInterestedMods = db_NSearch // Use the same glycan database as O-pair search
                .Where(b => b.ModificationType.Equals("N-linked glycosylation"))
                .Select(b => (b.ModificationType, b.IdWithMotif))
                .ToList();
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task_ModPair)}, new List<string> { raw }, new List<DbForTask>
            {
                db
            }, Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData")).Run();

            // Read ModSearch results
            var peptideResult_ModSearch = File.ReadLines(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData", "Task", "AllPSMs.psmtsv")).Skip(1).First().Split('\t');

            for (int i = 0; i < peptideResult_NSearch.Length; i++)
            {
                if (i == 13) // The full sequence of the peptide, there is little difference between the two searches
                {
                    Assert.That(peptideResult_NSearch[i].Equals("DAN[N-linked glycosylation:H5N2 on AnyWhere]NTQFQFTSR"));
                    Assert.That(peptideResult_ModSearch[i].Equals("DAN[N-linked glycosylation:H5N2 on Nxt]NTQFQFTSR"));
                }

                else
                {
                    Assert.That(peptideResult_NSearch[i], Is.EqualTo(peptideResult_ModSearch[i]), $"{peptideHearList[i]} don't match: expect {peptideResult_NSearch[i]} actual {peptideResult_ModSearch[i]}");
                }
            }
            Directory.Delete(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"), true);
        }

        [Test]
        public static void TestModPair_RegularMod()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSearchTaskconfigOGlycoTest_Run.toml"), MetaMorpheusTask.tomlConfig);

            DbForTask db = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P16150.fasta"), false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) }, new List<string> { spectraFile }, new List<DbForTask> { db }, outputFolder).Run();

            Directory.Delete(outputFolder, true);
        }

    }
}
