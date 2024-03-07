using EngineLayer;
using Nett;
using NUnit.Framework;
using Omics.Fragmentation;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using Omics.Fragmentation.Peptide;
using TaskLayer;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
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

        [Test]
        public static void CustomFragmentationManyTasks()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestCustomFragmentations");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");

            // setup tasks with custom ions
            var customIons = new List<ProductType>
                { ProductType.b, ProductType.c, ProductType.zDot };
            DissociationTypeCollection.ProductsFromDissociationType[DissociationType.Custom] = customIons;
            CalibrationTask calibrationTask = new CalibrationTask();
            calibrationTask.CommonParameters.GetType().GetProperty("DissociationType")?.SetValue(calibrationTask.CommonParameters, DissociationType.Custom);

            DissociationTypeCollection.ProductsFromDissociationType[DissociationType.Custom] = customIons;
            GptmdTask gptmdTask = new GptmdTask();
            gptmdTask.CommonParameters.GetType().GetProperty("DissociationType")?.SetValue(gptmdTask.CommonParameters, DissociationType.Custom);

            SearchTask searchTask = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBCZ.toml"),
                MetaMorpheusTask.tomlConfig);

            var taskCollection = new List<(string, MetaMorpheusTask)>
                { ("Calibration", calibrationTask), ("GPTMD", gptmdTask), ("Search", searchTask) };

            // assert custom ions are correct in common parameters for each task
            foreach (var task in taskCollection)
            {
                CollectionAssert.AreEquivalent(customIons, task.Item2.CommonParameters.CustomIons);
                Assert.That(task.Item2.CommonParameters.DissociationType, Is.EqualTo(DissociationType.Custom));
            }

            // run all tasks
            List<DbForTask> database = new List<DbForTask> { new DbForTask(myDatabase, false) };
            var engine = new EverythingRunnerEngine(taskCollection, new List<string> { myFile },
                database, outputFolder);
            engine.Run();

            // assert all file specific parameters contain the correct custom ions
            foreach (var task in taskCollection)
            {
                var fileSpecificCustomIons = task.Item2.FileSpecificParameters.First().Parameters.CustomIons;
                CollectionAssert.AreEqual(customIons, fileSpecificCustomIons);
            }

            // read generated toml settings and make sure custom fragmentations match are handled properly for each task
            var loadedCalibrationTask = Toml.ReadFile<CalibrationTask>(
                Path.Combine(outputFolder, "Task Settings", $"{taskCollection[0].Item1}Config.toml"),
                MetaMorpheusTask.tomlConfig);
            CollectionAssert.AreEquivalent(customIons, loadedCalibrationTask.CommonParameters.CustomIons);
            Assert.That(loadedCalibrationTask.CommonParameters.DissociationType, Is.EqualTo(DissociationType.Custom));
            var loadedGptmdTask = Toml.ReadFile<GptmdTask>(
                Path.Combine(outputFolder, "Task Settings", $"{taskCollection[1].Item1}Config.toml"),
                MetaMorpheusTask.tomlConfig);
            CollectionAssert.AreEquivalent(customIons, loadedGptmdTask.CommonParameters.CustomIons);
            Assert.That(loadedGptmdTask.CommonParameters.DissociationType, Is.EqualTo(DissociationType.Custom));
            var loadedSearchTask = Toml.ReadFile<SearchTask>(
                Path.Combine(outputFolder, "Task Settings", $"{taskCollection[2].Item1}Config.toml"),
                MetaMorpheusTask.tomlConfig);
            CollectionAssert.AreEquivalent(customIons, loadedSearchTask.CommonParameters.CustomIons);
            Assert.That(loadedSearchTask.CommonParameters.DissociationType, Is.EqualTo(DissociationType.Custom));

            // read gptmd and search results to ensure matched ions are correct
            var gptmdResults = PsmTsvReader.ReadTsv(Path.Combine(outputFolder, "GPTMD", "GPTMD_Candidates.psmtsv"), out List<string> warnings);
            Assert.That(!warnings.Any());
            var productIons = gptmdResults.SelectMany(p => p.MatchedIons.Select(m => m.NeutralTheoreticalProduct.ProductType))
                .Distinct();
            CollectionAssert.AreEquivalent(customIons, productIons);

            var searchResults = PsmTsvReader.ReadTsv(Path.Combine(outputFolder, "Search", "AllPSMs.psmtsv"), out warnings);
            Assert.That(!warnings.Any());
            productIons = searchResults.SelectMany(p => p.MatchedIons.Select(m => m.NeutralTheoreticalProduct.ProductType))
                .Distinct();
            CollectionAssert.AreEquivalent(customIons, productIons);

            // reload gptmd and search tasks with just calibrated files
            var newOutputFolder = Path.Combine(outputFolder, "ReSearch");
            var calibratedSpectraPath = Path.Combine(outputFolder, "Calibration", "SmallCalibratible_Yeast-calib.mzML");
            engine = new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)>() { ("GPTMD", loadedGptmdTask), ("Search", loadedSearchTask) },
                new List<string>() { calibratedSpectraPath },
                new List<DbForTask>() { new DbForTask(myDatabase, false) }, newOutputFolder);
            engine.Run();

            // read generated toml settings and make sure custom fragmentations match are handled properly for each task
            loadedGptmdTask = Toml.ReadFile<GptmdTask>(
                Path.Combine(newOutputFolder, "Task Settings", $"{taskCollection[1].Item1}Config.toml"),
                MetaMorpheusTask.tomlConfig);
            CollectionAssert.AreEquivalent(customIons, loadedGptmdTask.CommonParameters.CustomIons);
            Assert.That(loadedGptmdTask.CommonParameters.DissociationType, Is.EqualTo(DissociationType.Custom));
            loadedSearchTask = Toml.ReadFile<SearchTask>(
                Path.Combine(newOutputFolder, "Task Settings", $"{taskCollection[2].Item1}Config.toml"),
                MetaMorpheusTask.tomlConfig);
            CollectionAssert.AreEquivalent(customIons, loadedSearchTask.CommonParameters.CustomIons);
            Assert.That(loadedSearchTask.CommonParameters.DissociationType, Is.EqualTo(DissociationType.Custom));

            // read gptmd and search results to ensure matched ions are correct
            gptmdResults = PsmTsvReader.ReadTsv(Path.Combine(newOutputFolder, "GPTMD", "GPTMD_Candidates.psmtsv"), out warnings);
            Assert.That(!warnings.Any());
            productIons = gptmdResults.SelectMany(p => p.MatchedIons.Select(m => m.NeutralTheoreticalProduct.ProductType))
                .Distinct();
            CollectionAssert.AreEquivalent(customIons, productIons);

            searchResults = PsmTsvReader.ReadTsv(Path.Combine(newOutputFolder, "Search", "AllPSMs.psmtsv"), out warnings);
            Assert.That(!warnings.Any());
            productIons = searchResults.SelectMany(p => p.MatchedIons.Select(m => m.NeutralTheoreticalProduct.ProductType))
                .Distinct();
            CollectionAssert.AreEquivalent(customIons, productIons);

            // clean up
            Directory.Delete(outputFolder, true);

            // extra small tests to boost code coverage
            Assert.That(database.Count == 1);
            Assert.That(database.First().FileName, Is.EqualTo(Path.GetFileName(myDatabase)));

            loadedSearchTask.SearchParameters.CustomMdac = "Tacos";
            Assert.That(loadedSearchTask.SearchParameters.CustomMdac, Is.EqualTo("Tacos"));

            SingleTaskEventArgs args = new SingleTaskEventArgs("TaskTest");
            Assert.That(args.DisplayName, Is.EqualTo("TaskTest"));

            XmlForTaskListEventArgs xmlArgs =
                new XmlForTaskListEventArgs(new List<DbForTask> { new DbForTask(myDatabase, false) });
            Assert.That(xmlArgs.NewDatabases.Count == 1);
        }
    }
}