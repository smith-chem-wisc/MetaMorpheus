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
            string myFile1 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_A549_3_snip.mzML");
            string myFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_A549_3_snip_2.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_A549_3_snip.fasta");

            EverythingRunnerEngine engineToml = new(new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) }, new List<string> { myFile1, myFile2 }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engineToml.Run();

            string allResultsFile = Path.Combine(outputFolder, "allResults.txt");
            string[] allResults = File.ReadAllLines(allResultsFile);

            // TaGe_SA_A549_3_snip_2 is searched twice. First with two files being searched simultaneously, then with TaGe_SA_A549_3_snip_2 by itself
            // This allows us to compare the file specific results produced by in the two file search to the output
            // produced by searching the file by itself. The number of PSMs and Peptides observed should be the same
            // for both single-file and multi-file searches.
            // The number of protein groups will be different, because protein inference is performed once, using every peptide
            // identified across all files.
            int TaGe_SA_A549_3_snip_2ExpectedPsms = 215;
            int TaGe_SA_A549_3_snip_2ExpectedPeptides = 180;

            // The new PEP calculation method improves things, so all these numbers are increasing as of (7/17/24)
            Assert.AreEqual("All target PSMs with q-value = 0.01: 431", allResults[10]);
            Assert.AreEqual("All target peptides with q-value = 0.01: 180", allResults[11]);
            Assert.AreEqual("All target protein groups with q-value = 0.01 (1% FDR): 165", allResults[12]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - All target PSMs with q-value = 0.01: 216", allResults[14]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - All target peptides with q-value = 0.01: 170", allResults[15]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - All target PSMs with q-value = 0.01: 215", allResults[18]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - All target peptides with q-value = 0.01: 170", allResults[19]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target protein groups within 1 % FDR: 165", allResults[20]);

            string resultsFile = Path.Combine(outputFolder, @"postSearchAnalysisTaskTestOutput\results.txt");
            string[] results = File.ReadAllLines(resultsFile);
            Assert.AreEqual("All target PSMs with q-value = 0.01: 431", results[5]);
            Assert.AreEqual("All target peptides with q-value = 0.01: 180", results[6]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - All target PSMs with q-value = 0.01: 216", results[9]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - All target peptides with q-value = 0.01: 170", results[10]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target protein groups within 1 % FDR: 165", results[11]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - All target PSMs with q-value = 0.01: 215", results[13]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - All target peptides with q-value = 0.01: 170", results[14]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target protein groups within 1 % FDR: 165", results[15]);

            Directory.Delete(outputFolder, true);

            // Search TaGe_SA_A549_3_snip_2 by itself. The results from this should be identical to the file specific results above
            engineToml = new(new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) },
                new List<string> { myFile2 },
                new List<DbForTask> { new DbForTask(myDatabase, false) },
                outputFolder);
            engineToml.Run();

            string[] singleFileResults = File.ReadAllLines(resultsFile);
            Assert.AreEqual("All target PSMs with q-value = 0.01: " + TaGe_SA_A549_3_snip_2ExpectedPsms, singleFileResults[5]);
            Assert.AreEqual("All target peptides with q-value = 0.01: " + TaGe_SA_A549_3_snip_2ExpectedPeptides, singleFileResults[6]);
            Assert.AreEqual("All target protein groups with q-value = 0.01 (1% FDR): 165", singleFileResults[7]);

            //Second test that AllResults and Results display correct numbers of peptides and psms with PEP q-value filter on
            myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Task2-SearchTaskconfig.toml");
            searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            engineToml = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) }, new List<string> { myFile1, myFile2 }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engineToml.Run();

            allResultsFile = Path.Combine(outputFolder, "allResults.txt");
            allResults = File.ReadAllLines(allResultsFile);
            Assert.AreEqual("All target PSMs with pep q-value = 0.01: 528", allResults[10]);
            Assert.AreEqual("All target peptides with pep q-value = 0.01: 242", allResults[11]);
            Assert.AreEqual("All target protein groups with q-value = 0.01 (1% FDR): 138", allResults[12]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - All target PSMs with pep q-value = 0.01: 264", allResults[14]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - All target peptides with pep q-value = 0.01: 244", allResults[15]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target protein groups within 1 % FDR: 137", allResults[16]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - All target PSMs with pep q-value = 0.01: 264", allResults[18]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - All target peptides with pep q-value = 0.01: 308", allResults[19]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target protein groups within 1 % FDR: 138", allResults[20]);



            resultsFile = Path.Combine(outputFolder, @"postSearchAnalysisTaskTestOutput\results.txt");
            results = File.ReadAllLines(resultsFile);
            Assert.AreEqual("All target PSMs with pep q-value = 0.01: 528", results[5]);
            Assert.AreEqual("All target peptides with pep q-value = 0.01: 242", results[6]);
            Assert.AreEqual("All target protein groups with q-value = 0.01 (1% FDR): 138", results[7]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - All target PSMs with pep q-value = 0.01: 264", results[9]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - All target peptides with pep q-value = 0.01: 244", results[10]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target protein groups within 1 % FDR: 137", results[11]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - All target PSMs with pep q-value = 0.01: 264", results[13]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - All target peptides with pep q-value = 0.01: 308", results[14]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target protein groups within 1 % FDR: 138", results[15]);

            Directory.Delete(outputFolder, true);
        }
    }
}