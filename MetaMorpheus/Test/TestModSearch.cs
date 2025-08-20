using Easy.Common.Extensions;
using EngineLayer;
using EngineLayer.GlycoSearch;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics;
using Nett;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using Omics;
using Omics.Digestion;
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
using System.Windows.Documents;
using System.Windows.Media.Animation;
using TaskLayer;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;

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
            List<string> diffList = new List<string>();
            for (int i = 0; i < 37; i++)
            {
                if (i == 28 || i == 31 || i == 32)
                {
                    continue;
                }

                if (peptideResult_ModSearch[i] != peptideResult_NSearch[i])
                {
                    diffList.Add(peptideHearList[i]);
                }

                //Assert.That(peptideResult_NSearch[i], Is.EqualTo(peptideResult_ModSearch[i]), $"{peptideHearList[i]} don't match: expect {peptideResult_NSearch[i]} actual {peptideResult_ModSearch[i]}");
            }
            Directory.Delete(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"), true);
        }

        [Test]
        public static void TestModPair_GlobalMod()
        {
            //This test will compare the search performance on GPTMD and ModPair while using the global mod list.
            string outputFolder = "E:\\ModPair\\RegularCompare\\ModPair_EtHCD";

            //Create Search Task
            string fasta = "E:\\ModPair\\RegularCompare\\uniprotkb_AND_reviewed_true_AND_model_o_2025_08_18.fasta";
            DbForTask db = new DbForTask(fasta, false);
            List<(string, string)> interestMods = new List<(string, string)>()
            {
                ("Common Biological","Phosphorylation on S"), ("Common Biological","Phosphorylation on T")
            };

            ////Create the spectra file
            //string xml = "E:\\ModPair\\RegularCompare\\uniprotkb_AND_reviewed_true_AND_model_o_2025_08_18.xml";
            //var proteins = ProteinDbLoader.LoadProteinXML(xml, true, DecoyType.None, GlobalVariables.AllModsKnown, false, null,
            //    out _);

            //var peptides = proteins
            //    .SelectMany(p => p.Digest((IDigestionParams) new DigestionParams(), new List<Modification>(), new List<Modification>()))
            //    .ToList();
            //var modPeptideList = peptides.Where(p => p.AllModsOneIsNterminus.Count >= 1)
            //    .ToList();
            //string[] validModNames = ["Phospho"];
            //var filteredModifiedPeptideList = modPeptideList.Where(p =>
            //    p.AllModsOneIsNterminus.All(mod =>
            //        validModNames.Any(validName => mod.Value.IdWithMotif.Contains(validName, StringComparison.InvariantCultureIgnoreCase)))).Take(1000).ToList();

            //string rawFile = "fakeRaw_EtHCD.mzML";
            //MsDataFile fakeRawFile = new TestDataFile(filteredModifiedPeptideList);
            //Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(fakeRawFile, rawFile, false);
            //string peptidePath = "E:\\ModPair\\RegularCompare\\peptidelist";
            //File.WriteAllLines(peptidePath, filteredModifiedPeptideList.Select(p=>p.FullSequence).ToList());


            // First, we run the GPTMD search
            //var engine = new EverythingRunnerEngine(
            //    new List<(string, MetaMorpheusTask)> { ("task1", task1) },
            //    new List<string>() { rawFile },
            //    new List<DbForTask> { new DbForTask(proteinDB, false) },
            //    Environment.CurrentDirectory);


            // Second, we run the ModPair search


            //Load the toml file and modify the parameters
            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSearchTaskconfigModPairTest_Run.toml"), MetaMorpheusTask.tomlConfig);

            glycoSearchTask._glycoSearchParameters.GlycoSearchType = GlycoSearchType.ModSearch;
            glycoSearchTask._glycoSearchParameters.MaximumOGlycanAllowed = 2;
            glycoSearchTask._glycoSearchParameters.ListOfInterestedMods = interestMods;
            string spectraFile = "E:\\ModPair\\RegularCompare\\fakeRaw_EtHCD.mzML";
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) }, new List<string> { spectraFile }, new List<DbForTask> { db }, outputFolder).Run();
            int iii = 0;
        }

        [Test]
        public static void ComparePsm()
        {
            string anserPsm = "E:\\ModPair\\RegularCompare\\peptidelist.txt";
            string metaPsm = "E:\\ModPair\\RegularCompare\\MetaSearch_EtHCD\\Task1-Search\\AllPSMs.psmtsv";
            string modPairPsm = "E:\\ModPair\\RegularCompare\\ModPair_EtHCD\\Task\\AllPSMs.psmtsv";
            string GptmdPsm = "E:\\ModPair\\RegularCompare\\GPTMD_EtHCD\\Task2-Search\\AllPSMs.psmtsv";

            HashSet<string> answerSet = new HashSet<string>();
            foreach (var line in File.ReadLines(anserPsm))
            {
                string newline = line;
                if (newline.Contains("[UniProt:Phosphoserine on S]"))
                {
                    newline = newline.Replace("[UniProt:Phosphoserine on S]", "[Common Biological:Phosphorylation on S]");
                }
                if (newline.Contains("[UniProt:Phosphothreonine on T]"))
                {
                    newline = newline.Replace("[UniProt:Phosphothreonine on T]", "[Common Biological:Phosphorylation on T]");
                }
                answerSet.Add(newline);
            }
            HashSet<string> metaSet = new HashSet<string>();
            foreach (var line in File.ReadLines(metaPsm).Skip(0))
            {
                var lines = line.Split('\t');
                if (lines[13] == lines[14])
                {
                    continue; 
                }
                if (line.Split('\t')[33].Contains('Y') || line.Split('\t')[32].Contains('Y'))
                {
                    continue;
                }

                string newline = line.Split('\t')[14];
                if (newline.Contains("[UniProt:Phosphoserine on S]"))
                {
                    newline = newline.Replace("[UniProt:Phosphoserine on S]", "[Common Biological:Phosphorylation on S]");
                }
                if (newline.Contains("[UniProt:Phosphothreonine on T]"))
                {
                    newline = newline.Replace("[UniProt:Phosphothreonine on T]", "[Common Biological:Phosphorylation on T]");
                }
                metaSet.Add(newline);
            }

            HashSet<string> modPairSet = new HashSet<string>();
            foreach (var line in File.ReadLines(modPairPsm))
            {
                var lines = line.Split('\t');
                if (lines[11] == lines[13])
                {
                    continue;
                }
                if (lines[24] == "D")
                {
                    continue;
                }

                string full = line.Split('\t')[13];
                modPairSet.Add(full);
            }
            HashSet<string> gptmdSet = new HashSet<string>();
            foreach (var line in File.ReadLines(GptmdPsm).Skip(0))
            {
                var lines = line.Split('\t');
                if (lines[13] == lines[14])
                {
                    continue;
                }
                if (lines[33].Contains('Y') || lines[32].Contains('Y'))
                {
                    continue;
                }
                string full = line.Split('\t')[14];
                gptmdSet.Add(full);
            }

            // Find strings in answerSet but not in metaSet
            var notInMeta = answerSet.Except(metaSet).ToList();
            // Find strings in answerSet but not in modPairSet
            var notInModPair = answerSet.Except(modPairSet).ToList();
            // Find strings in answerSet but not in gptmdSet
            var notInGptmd = answerSet.Except(gptmdSet).ToList();

            int iiii = 0;
        }

        [Test]
        public static void TestSub()
        {
            string outputFolder = "E:\\ModPair\\zh\\test1";

            string raw = "E:\\ModPair\\zh\\test1\\01-14-25_bu-DDA_Yeast_SP3_sv95.mzML";
            string fasta = "E:\\ModPair\\zh\\test1\\modifiedYeastFasta100_noAmb_noMiss_new_100only.fasta";
            DbForTask db = new DbForTask(fasta, false);
            string GPTMDToml = "E:\\ModPair\\zh\\test1\\Task1-GPTMDTaskconfig.toml";
            string search = "E:\\ModPair\\zh\\test1\\Task2-SearchTaskconfig.toml";
            string modPair = "E:\\ModPair\\zh\\test1\\GlycoSearchTaskconfigModPairTest_Run.toml";

            var task1 = Toml.ReadFile<GptmdTask>(GPTMDToml, MetaMorpheusTask.tomlConfig);

            List<(string, string)> interestMods = task1.GptmdParameters.ListOfModsGptmd;
            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(modPair, MetaMorpheusTask.tomlConfig);
            glycoSearchTask._glycoSearchParameters.GlycoSearchType = GlycoSearchType.ModSearch;
            glycoSearchTask._glycoSearchParameters.MaximumOGlycanAllowed = 1;
            glycoSearchTask._glycoSearchParameters.ListOfInterestedMods = interestMods;
            glycoSearchTask._glycoSearchParameters.OxoniumIonFilt = false;

            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) }, new List<string> { raw }, new List<DbForTask> { db }, outputFolder).Run();
            int iii = 0;

        }



    }
}
