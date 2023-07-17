using System.Collections.Generic;
using System.IO;
using Nett;
using NUnit.Framework;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class PostSearchAnalysisTaskTests
    {
        [Test]
        public static void AllResultsAndResultsTxtTests()
        {
            //First test that AllResults and Results display correct numbers of peptides and psms with q-value filter on
            string myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Task1-SearchTaskconfig.toml");
            SearchTask searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PostSearchAnalysisTaskTest");
            string myFile1 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\miniA549.mzML");
            string myFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\miniA549_2.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\miniA549_2.fasta");

            EverythingRunnerEngine engineToml = new(new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) }, new List<string> { myFile1, myFile2 }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engineToml.Run();

            string allResultsFile = Path.Combine(outputFolder, "allResults.txt");
            string[] allResults = File.ReadAllLines(allResultsFile);

            // b6_2 is searched twice. First with two files being searched simultaneously, then with b6_2 by itself
            // This allows us to compare the file specific results produced by in the two file search to the output
            // produced by searching the file by itself. The number of PSMs and Peptides observed should be the same
            // for both single-file and multi-file searches.
            // The number of protein groups will be different, because protein inference is performed once, using every peptide
            // identified across all files.

            Assert.AreEqual("All target PSMs with q-value = 0.01: 228", allResults[12]);
            Assert.AreEqual("All target peptides with q-value = 0.01 : 178", allResults[13]);
            Assert.AreEqual("All target protein groups with q-value = 0.01 (1% FDR): 170", allResults[14]);
            Assert.AreEqual("miniA549_2 target PSMs with q-value = 0.01: 114", allResults[18]);
            Assert.AreEqual("Target protein groups within 1 % FDR in miniA549_2: 95", allResults[24]);
            Assert.AreEqual("miniA549_2 Target peptides with q-value = 0.01 : 0", allResults[26]);
            
            Assert.AreEqual("miniA549 target PSMs with q-value = 0.01: 111", allResults[22]);
            Assert.AreEqual("miniA549 Target peptides with q-value = 0.01 : 103", allResults[28]);
            Assert.AreEqual("Target protein groups within 1 % FDR in miniA549: 97", allResults[25]);

            string resultsFile = Path.Combine(outputFolder, @"postSearchAnalysisTaskTestOutput\results.txt");
            string[] results = File.ReadAllLines(resultsFile);
            Assert.AreEqual("All target PSMs with q-value = 0.01: 228", results[7]);
            Assert.AreEqual("All target peptides with q-value = 0.01 : 178", results[8]);
            Assert.AreEqual("All target protein groups with q-value = 0.01 (1% FDR): 170", results[9]);
            Assert.AreEqual("miniA549_2 target PSMs with q-value = 0.01: 114", results[13]);
            Assert.AreEqual("Target protein groups within 1 % FDR in miniA549_2: 95", results[19]);
            Assert.AreEqual("miniA549_2 Target peptides with q-value = 0.01 : 0", results[21]);
            Assert.AreEqual("Target protein groups within 1 % FDR in miniA549: 97", results[20]);

            Assert.AreEqual("miniA549 target PSMs with q-value = 0.01: 111", results[17]);
            Assert.AreEqual("miniA549 Target peptides with q-value = 0.01 : 103", results[23]);
            

            Directory.Delete(outputFolder, true);

            // Search slice_b6_2 by itself. The results from this should be identical to the file specific results above
            engineToml = new(new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) },
                new List<string> { myFile2 },
                new List<DbForTask> { new DbForTask(myDatabase, false) },
                outputFolder);
            engineToml.Run();

            string[] singleFileResults = File.ReadAllLines(resultsFile);
            Assert.AreEqual("All target PSMs with q-value = 0.01: 114", singleFileResults[7]); 
            Assert.AreEqual("All target peptides with q-value = 0.01 : 0", singleFileResults[8]); 
            Assert.AreEqual("All target protein groups with q-value = 0.01 (1% FDR): 92", singleFileResults[9]); 

            //Second test that AllResults and Results display correct numbers of peptides and psms with PEP q-value filter on
            myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Task2-SearchTaskconfig.toml");
            searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            engineToml = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) }, new List<string> { myFile1, myFile2 }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engineToml.Run();

            allResultsFile = Path.Combine(outputFolder, "allResults.txt");
            allResults = File.ReadAllLines(allResultsFile);
            Assert.AreEqual("All target PSMs with pep q-value = 0.01: 278", allResults[12]);
            Assert.AreEqual("All target peptides with pep q-value = 0.01 : 210", allResults[13]);
            Assert.AreEqual("All target protein groups with q-value = 0.01 (1% FDR): 170", allResults[14]);
            Assert.AreEqual("miniA549_2 target PSMs with pep q-value = 0.01: 152", allResults[18]);
            Assert.AreEqual("miniA549 target PSMs with pep q-value = 0.01: 125", allResults[22]);
            Assert.AreEqual("Target protein groups within 1 % FDR in miniA549_2: 95", allResults[24]);
            Assert.AreEqual("Target protein groups within 1 % FDR in miniA549: 97", allResults[25]);
            Assert.AreEqual("miniA549_2 Target peptides with pep q-value = 0.01 : 125", allResults[26]);
            Assert.AreEqual("miniA549 Target peptides with pep q-value = 0.01 : 112", allResults[28]);


            resultsFile = Path.Combine(outputFolder, @"postSearchAnalysisTaskTestOutput\results.txt");
            results = File.ReadAllLines(resultsFile);
            Assert.AreEqual("All target PSMs with pep q-value = 0.01: 278", results[7]);
            Assert.AreEqual("All target peptides with pep q-value = 0.01 : 210", results[8]);
            Assert.AreEqual("All target protein groups with q-value = 0.01 (1% FDR): 170", results[9]);
            Assert.AreEqual("miniA549_2 target PSMs with pep q-value = 0.01: 152", results[13]);
            Assert.AreEqual("miniA549 target PSMs with pep q-value = 0.01: 125", results[17]);
            Assert.AreEqual("Target protein groups within 1 % FDR in miniA549_2: 95", results[19]);
            Assert.AreEqual("Target protein groups within 1 % FDR in miniA549: 97", results[20]);
            Assert.AreEqual("miniA549_2 Target peptides with pep q-value = 0.01 : 125", results[21]);
            Assert.AreEqual("miniA549 Target peptides with pep q-value = 0.01 : 112", results[23]);

            Directory.Delete(outputFolder, true);
        }
    }
}
