﻿using EngineLayer;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class TestToml
    {
        [Test]
        public static void TestTomlFunction()
        {
            SearchTask searchTask = new SearchTask
            {
                CommonParameters = new CommonParameters(
                    productMassTolerance: new PpmTolerance(666),
                    listOfModsFixed: new List<(string, string)> { ("a", "b"), ("c", "d") },
                    listOfModsVariable: new List<(string, string)> { ("e", "f"), ("g", "h") }),
            };
            Toml.WriteFile(searchTask, "SearchTask.toml", MetaMorpheusTask.tomlConfig);
            var searchTaskLoaded = Toml.ReadFile<SearchTask>("SearchTask.toml", MetaMorpheusTask.tomlConfig);

            Assert.AreEqual(searchTask.CommonParameters.DeconvolutionMassTolerance.ToString(), searchTaskLoaded.CommonParameters.DeconvolutionMassTolerance.ToString());
            Assert.AreEqual(searchTask.CommonParameters.ProductMassTolerance.ToString(), searchTaskLoaded.CommonParameters.ProductMassTolerance.ToString());
            Assert.AreEqual(searchTask.CommonParameters.PrecursorMassTolerance.ToString(), searchTaskLoaded.CommonParameters.PrecursorMassTolerance.ToString());

            Assert.AreEqual(searchTask.CommonParameters.ListOfModsFixed.Count(), searchTaskLoaded.CommonParameters.ListOfModsFixed.Count());
            Assert.AreEqual(searchTask.CommonParameters.ListOfModsFixed.First().Item1, searchTaskLoaded.CommonParameters.ListOfModsFixed.First().Item1);
            Assert.AreEqual(searchTask.CommonParameters.ListOfModsFixed.First().Item2, searchTaskLoaded.CommonParameters.ListOfModsFixed.First().Item2);

            Assert.AreEqual(searchTask.CommonParameters.ListOfModsVariable.Count(), searchTaskLoaded.CommonParameters.ListOfModsVariable.Count());

            Assert.AreEqual(searchTask.SearchParameters.MassDiffAcceptorType, searchTaskLoaded.SearchParameters.MassDiffAcceptorType);
            Assert.AreEqual(searchTask.SearchParameters.CustomMdac, searchTaskLoaded.SearchParameters.CustomMdac);

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestConsistency");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");

            var engine = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Search", searchTask) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engine.Run();
            var engineToml = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("SearchTOML", searchTaskLoaded) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engineToml.Run();

            var results = File.ReadAllLines(Path.Combine(outputFolder, @"Search\AllPSMs.psmtsv"));
            var resultsToml = File.ReadAllLines(Path.Combine(outputFolder, @"SearchTOML\AllPSMs.psmtsv"));
            Assert.That(results.SequenceEqual(resultsToml));

            CalibrationTask calibrationTask = new CalibrationTask();
            Toml.WriteFile(calibrationTask, "CalibrationTask.toml", MetaMorpheusTask.tomlConfig);
            var calibrationTaskLoaded = Toml.ReadFile<CalibrationTask>("CalibrationTask.toml", MetaMorpheusTask.tomlConfig);

            GptmdTask gptmdTask = new GptmdTask();
            Toml.WriteFile(gptmdTask, "GptmdTask.toml", MetaMorpheusTask.tomlConfig);
            var gptmdTaskLoaded = Toml.ReadFile<GptmdTask>("GptmdTask.toml", MetaMorpheusTask.tomlConfig);

            var gptmdEngine = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("GPTMD", gptmdTask) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            gptmdEngine.Run();
            var gptmdEngineToml = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("GPTMDTOML", gptmdTaskLoaded) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            gptmdEngineToml.Run();

            var gptmdResults = File.ReadAllLines(Path.Combine(outputFolder, @"GPTMD\GPTMD_Candidates.psmtsv"));
            var gptmdResultsToml = File.ReadAllLines(Path.Combine(outputFolder, @"GPTMDTOML\GPTMD_Candidates.psmtsv"));

            Assert.That(gptmdResults.SequenceEqual(gptmdResultsToml));

            XLSearchTask xLSearchTask = new XLSearchTask();
            Toml.WriteFile(xLSearchTask, "XLSearchTask.toml", MetaMorpheusTask.tomlConfig);
            var xLSearchTaskLoaded = Toml.ReadFile<XLSearchTask>("XLSearchTask.toml", MetaMorpheusTask.tomlConfig);

            string myFileXl = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA_DSSO_ETchD6010.mgf");
            string myDatabaseXl = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA.fasta");

            var xlEngine = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("XLSearch", xLSearchTask) }, new List<string> { myFileXl }, new List<DbForTask> { new DbForTask(myDatabaseXl, false) }, outputFolder);
            xlEngine.Run();
            var xlEngineToml = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("XLSearchTOML", xLSearchTaskLoaded) }, new List<string> { myFileXl }, new List<DbForTask> { new DbForTask(myDatabaseXl, false) }, outputFolder);
            xlEngineToml.Run();

            var xlResults = File.ReadAllLines(Path.Combine(outputFolder, @"XLSearch\XL_Intralinks.tsv"));
            var xlResultsToml = File.ReadAllLines(Path.Combine(outputFolder, @"XLSearchTOML\XL_Intralinks.tsv"));

            Assert.That(xlResults.SequenceEqual(xlResultsToml));
            Directory.Delete(outputFolder, true);
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GptmdTask.toml"));
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"XLSearchTask.toml"));
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"SearchTask.toml"));
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"CalibrationTask.toml"));
        }

        [Test]
        public static void TestTomlForSpecficFiles()
        {
            var fileSpecificToml = Toml.ReadFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "testFileSpecfic.toml"), MetaMorpheusTask.tomlConfig);
            var tomlSettingsList = fileSpecificToml.ToDictionary(p => p.Key);
            Assert.AreEqual(tomlSettingsList["Protease"].Value.Get<string>(), "Asp-N");
            Assert.IsFalse(tomlSettingsList.ContainsKey("maxMissedCleavages"));
            Assert.IsFalse(tomlSettingsList.ContainsKey("InitiatorMethionineBehavior"));

            FileSpecificParameters f = new FileSpecificParameters(fileSpecificToml);

            Assert.AreEqual("Asp-N", f.Protease.Name);
            Assert.IsNull(f.MaxMissedCleavages);

            CommonParameters c = MetaMorpheusTask.SetAllFileSpecificCommonParams(new CommonParameters(), f);

            Assert.AreEqual("Asp-N", c.DigestionParams.Protease.Name);
            Assert.AreEqual(2, c.DigestionParams.MaxMissedCleavages);
        }
    }
}