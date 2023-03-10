using EngineLayer;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using Proteomics.Fragmentation;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class CustomFragmentationTest
    {
        [Test]
        public static void MultipleCustomFragmentations()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestCustomFragmentations");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");

            // create 3 search tasks with different custom fragmentation ions
            var task1 = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBY.toml"), MetaMorpheusTask.tomlConfig);
            var task2 = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customCZ.toml"), MetaMorpheusTask.tomlConfig);
            var task3 = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBCZ.toml"), MetaMorpheusTask.tomlConfig);

            // run all tasks
            DbForTask db = new DbForTask(myDatabase, false);
            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("TestSearchBY", task1), ("TestSearchCZ", task2), ("TestSearchBCZ", task3) };
            var engine = new EverythingRunnerEngine(taskList, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engine.Run();

            // read generated toml settings and make sure custom fragmentations match are handled properly for each task
            string outputPath = Path.Combine(outputFolder, @"Task Settings");

            var task1Settings = Toml.ReadFile<SearchTask>(Path.Combine(outputPath, @"TestSearchBYConfig.toml"), MetaMorpheusTask.tomlConfig);
            var task2Settings = Toml.ReadFile<SearchTask>(Path.Combine(outputPath, @"TestSearchCZConfig.toml"), MetaMorpheusTask.tomlConfig);
            var task3Settings = Toml.ReadFile<SearchTask>(Path.Combine(outputPath, @"TestSearchBCZConfig.toml"), MetaMorpheusTask.tomlConfig);

            var task1CustomIons = task1Settings.CommonParameters.CustomIons;
            var task2CustomIons = task2Settings.CommonParameters.CustomIons;
            var task3CustomIons = task3Settings.CommonParameters.CustomIons;

            Assert.That(task1CustomIons.Contains(ProductType.b) && task1CustomIons.Contains(ProductType.y) && task1CustomIons.Count == 2);
            Assert.That(task2CustomIons.Contains(ProductType.c) && task2CustomIons.Contains(ProductType.zDot) && task2CustomIons.Count == 2);
            Assert.That(task3CustomIons.Contains(ProductType.c) && task3CustomIons.Contains(ProductType.zDot) && task3CustomIons.Contains(ProductType.b) && task3CustomIons.Count == 3);

            Directory.Delete(outputFolder, true);
        }       
    }
}
