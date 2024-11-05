using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using EngineLayer.FdrAnalysis;
using Microsoft.VisualStudio.TestPlatform.ObjectModel;
using Nett;
using NUnit.Framework;
using TaskLayer;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

namespace Test
{
    /// <summary>
    /// Uses test cases found in EverythingRunnerEngineTestCase.cs
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class PostSearchAnalysisTaskTestsPsmCountOverride
    {

        [Test]
        [NonParallelizable]
        public static void AllResultsAndResultsTxtContainsCorrectValues_QValue_BottomUp()
        {
            //override to be only used for unit tests in non-parallelizable format
            //must set to false at the end of this method
            var type = typeof(FdrAnalysisEngine);
            var property = type.GetProperty("QvalueThresholdOverride");
            property.SetValue(null, true);

            var ResultDirectory =
            Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PsmCountOverrideTest");

            // Directory GlobalSetup
            if (Directory.Exists(ResultDirectory))
                Directory.Delete(ResultDirectory, true);

            if (!Directory.Exists(ResultDirectory))
                Directory.CreateDirectory(ResultDirectory);

            string myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\Task1-SearchTaskconfig.toml");
            SearchTask searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            string myFile1 = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\TaGe_SA_A549_3_snip.mzML");
            string myFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\TaGe_SA_A549_3_snip_2.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\TaGe_SA_A549_3_snip.fasta");
            searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            searchTaskLoaded.SearchParameters.WriteIndividualFiles = true;
            searchTaskLoaded.SearchParameters.WriteMzId = true;
            searchTaskLoaded.SearchParameters.WritePepXml = true;

            List<(string, MetaMorpheusTask)> taskList = new() { ("searchTask", searchTaskLoaded) };
            List<DbForTask> DatabaseList = new() { new DbForTask(myDatabase, false) };
            List<string> DataFileList = new() { myFile1, myFile2 };
            string OutputDirectory = Path.Combine(ResultDirectory, "BottomUpQValue");

            var e = new EverythingRunnerEngine(taskList, DataFileList, DatabaseList, OutputDirectory);

            e.Run();

            string allResultsFile = Path.Combine(OutputDirectory, "allResults.txt");
            string[] allResults = File.ReadAllLines(allResultsFile);

            // The new PEP calculation method improves things, so all these numbers are increasing as of (7/17/24)
            // There is a discrepancy between the number of All target peptides and individual file target peptides, 
            // presumably due to the way protein inference is performed.
            Assert.AreEqual("All target PSMs with q-value <= 0.01: 428", allResults[10]);
            Assert.AreEqual("All target peptides with q-value <= 0.01: 174", allResults[11]);
            Assert.AreEqual("All target protein groups with q-value <= 0.01 (1% FDR): 165", allResults[12]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target PSMs with q-value <= 0.01: 214", allResults[14]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target peptides with q-value <= 0.01: 174", allResults[15]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target protein groups within 1 % FDR: 165", allResults[16]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target PSMs with q-value <= 0.01: 214", allResults[18]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target peptides with q-value <= 0.01: 174", allResults[19]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target protein groups within 1 % FDR: 165", allResults[20]);


            string resultsFile = Path.Combine(OutputDirectory, @"searchTask\results.txt");
            string[] results = File.ReadAllLines(resultsFile);

            Assert.AreEqual("All target PSMs with q-value <= 0.01: 428", results[5]);
            Assert.AreEqual("All target peptides with q-value <= 0.01: 174", results[6]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target PSMs with q-value <= 0.01: 214", results[9]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target peptides with q-value <= 0.01: 174", results[10]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target protein groups within 1 % FDR: 165", results[11]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target PSMs with q-value <= 0.01: 214", results[13]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target peptides with q-value <= 0.01: 174", results[14]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target protein groups within 1 % FDR: 165", results[15]);

            property.SetValue(null, false);

            // Directory GlobalSetup
            if (Directory.Exists(ResultDirectory))
                Directory.Delete(ResultDirectory, true);
        }

        [Test]
        [NonParallelizable]
        public static void AllResultsAndResultsTxtContainsCorrectValues_PepQValue_BottomUp()
        {
            //override to be only used for unit tests in non-parallelizable format
            //must set to false at the end of this method
            var type = typeof(FdrAnalysisEngine);
            var property = type.GetProperty("QvalueThresholdOverride");
            property.SetValue(null, true);

            var ResultDirectory =
            Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PsmCountOverrideTest");

            // Directory GlobalSetup
            if (Directory.Exists(ResultDirectory))
                Directory.Delete(ResultDirectory, true);

            if (!Directory.Exists(ResultDirectory))
                Directory.CreateDirectory(ResultDirectory);

            string myFile1 = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\TaGe_SA_A549_3_snip.mzML");
            string myFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\TaGe_SA_A549_3_snip_2.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\TaGe_SA_A549_3_snip.fasta");

            var myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\Task2-SearchTaskconfig.toml");
            var searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            searchTaskLoaded.CommonParameters.QValueCutoffForPepCalculation = 0.01;
            searchTaskLoaded.SearchParameters.WriteIndividualFiles = true;
            searchTaskLoaded.SearchParameters.WriteMzId = true;
            searchTaskLoaded.SearchParameters.WritePepXml = true;

            List<(string, MetaMorpheusTask)> taskList = new() { ("searchTask", searchTaskLoaded) };
            List<DbForTask> DatabaseList = new() { new DbForTask(myDatabase, false) };
            List<string> DataFileList = new() { myFile1, myFile2 };
            string OutputDirectory = Path.Combine(ResultDirectory, "BottomUpPepQValue");

            var e = new EverythingRunnerEngine(taskList, DataFileList, DatabaseList, OutputDirectory);

            e.Run();

            string allResultsFile = Path.Combine(OutputDirectory, "allResults.txt");
            string[] allResults = File.ReadAllLines(allResultsFile);

            // The new PEP calculation method improves things, so all these numbers are increasing as of (7/17/24)
            // There is a discrepancy between the number of All target peptides and individual file target peptides, 
            // presumably due to the way protein inference is performed.
            Assert.AreEqual("All target PSMs with pep q-value <= 0.01: 382", allResults[10]);
            Assert.AreEqual("All target peptides with pep q-value <= 0.01: 153", allResults[11]);
            Assert.AreEqual("All target protein groups with q-value <= 0.01 (1% FDR): 140", allResults[12]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target PSMs with pep q-value <= 0.01: 190", allResults[14]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target peptides with pep q-value <= 0.01: 153", allResults[15]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target protein groups within 1 % FDR: 140", allResults[16]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target PSMs with pep q-value <= 0.01: 190", allResults[18]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target peptides with pep q-value <= 0.01: 153", allResults[19]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target protein groups within 1 % FDR: 140", allResults[20]);


            string resultsFile = Path.Combine(OutputDirectory, @"searchTask\results.txt");
            string[] results = File.ReadAllLines(resultsFile);

            Assert.AreEqual("All target PSMs with pep q-value <= 0.01: 382", results[5]);
            Assert.AreEqual("All target peptides with pep q-value <= 0.01: 153", results[6]);
            Assert.AreEqual("All target protein groups with q-value <= 0.01 (1% FDR): 140", results[7]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target PSMs with pep q-value <= 0.01: 190", results[9]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target peptides with pep q-value <= 0.01: 153", results[10]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target protein groups within 1 % FDR: 140", results[11]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target PSMs with pep q-value <= 0.01: 190", results[13]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target peptides with pep q-value <= 0.01: 153", results[14]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target protein groups within 1 % FDR: 140", results[15]);

            property.SetValue(null, false);

            // Directory GlobalSetup
            if (Directory.Exists(ResultDirectory))
                Directory.Delete(ResultDirectory, true);
        }

    }
}