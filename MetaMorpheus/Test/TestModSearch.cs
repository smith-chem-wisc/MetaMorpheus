using Easy.Common.Extensions;
using Nett;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Animation;
using EngineLayer.GlycoSearch;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public class TestModSearch
    {
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
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSearchTaskconfigOGlycoTest_Run.toml"), MetaMorpheusTask.tomlConfig);

            DbForTask db = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P16150.fasta"), false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) }, new List<string> { spectraFile }, new List<DbForTask> { db }, outputFolder).Run();

            Directory.Delete(outputFolder, true);
        }


        [Test]
        public static void TestModPair_N()
        {
            // set up output directory, MS data, and database
            Directory.CreateDirectory(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"));
                DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\Q9C0Y4.fasta"), false);
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\yeast_glycan_25170.mgf");

            // first, we run the NPair search (O-pair for Nglyco)
            var task_NPair = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\NGlycanSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task_NPair) }, new List<string> { raw }, new List<DbForTask>
            {
                db
            }, Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData")).Run();
            var peptideHearList = File
                .ReadLines(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData", "Task", "AllPSMs.psmtsv"))
                .First().Split('\t');
            var peptideResult_NSearch = File.ReadLines(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData", "Task", "AllPSMs.psmtsv"))
                .Skip(1).First().Split('\t');
            Directory.Delete(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"), true);



            // Now run the ModPair search
            Directory.CreateDirectory(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"));
            var task_ModPair = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\NGlycanSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            task_ModPair._glycoSearchParameters.GlycoSearchType = GlycoSearchType.ModSearch;
            task_ModPair._glycoSearchParameters.MaximumOGlycanAllowed = 1;
            
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task_ModPair)}, new List<string> { raw }, new List<DbForTask>
            {
                db
            }, Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData")).Run();
            // Read ModSearch results
            var peptideResult_ModSearch = File.ReadLines(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData", "Task", "AllPSMs.psmtsv")).Skip(1).First().Split('\t');
            for (int i = 0; i < peptideResult_NSearch.Length; i++)
            {
                if (i == 13)
                {
                    var fullSeq_NSearch = peptideResult_NSearch[i];
                    var fullSeq_ModSearch = peptideResult_ModSearch[i];
                    Assert.That(fullSeq_NSearch.Equals("DAN[N-linked glycosylation:H5N2 on AnyWhere]NTQFQFTSR"));
                    Assert.That(fullSeq_ModSearch.Equals("DAN[N-linked glycosylation:H5N2 on Nxt]NTQFQFTSR"));
                }

                else
                {
                    Assert.That(peptideResult_NSearch[i], Is.EqualTo(peptideResult_ModSearch[i]), $"{peptideHearList[i]} don't match");
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
