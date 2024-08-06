using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using Easy.Common.Extensions;
using EngineLayer;
using Nett;
using NUnit.Framework; using Assert = NUnit.Framework.Legacy.ClassicAssert;
using TaskLayer;

namespace Test
{
    public class SearchTaskTestCase
    {
        
        internal string Identifier { get; set; }
        internal string OutputDirectory { get; set; }
        internal int SpectraFileCount { get; set; }
        internal bool IsTopDown { get; init; }
        internal bool WriteIndividualResults { get; init; }
        internal bool WritePepXml { get; init; }
        internal bool WritePercolatorFiles { get; init; }

        #region TestCases
        internal static string ResultDirectory => Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PostSearchAnalysisTaskTest");

        internal static SearchTaskTestCase TestCaseLocator(string identifier)
        {
            return identifier switch
            {
                "BottomUpQValue" => BottomUpQValue,
                "BottomUpQValueSingle" => BottomUpQValueSingle,
                "BottomUpPepQValue" => BottomUpPepQValue,
                _ => throw new ArgumentOutOfRangeException(nameof(identifier), identifier, null)
            };
        }

        internal static string[] AllTestCaseIdentifiers => new[]
        {
            "BottomUpQValue",
            "BottomUpQValueSingle",
            "BottomUpPepQValue"
        };

        private static SearchTaskTestCase _bottomUpQValue;
        public static SearchTaskTestCase BottomUpQValue
        {

            get
            {
                if (!_bottomUpQValue.IsDefault()) return _bottomUpQValue;

                string outputFolder = Path.Combine(ResultDirectory, @"BottomUpQValue");
                if (Directory.Exists(outputFolder))
                    Directory.Delete(outputFolder, true);
                string myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                    @"TestData\Task1-SearchTaskconfig.toml");
                SearchTask searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
                string myFile1 = Path.Combine(TestContext.CurrentContext.TestDirectory,
                    @"TestData\TaGe_SA_A549_3_snip.mzML");
                string myFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory,
                    @"TestData\TaGe_SA_A549_3_snip_2.mzML");
                string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory,
                    @"TestData\TaGe_SA_A549_3_snip.fasta");
                EverythingRunnerEngine engineToml =
                    new(new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) },
                        new List<string> { myFile1, myFile2 }, new List<DbForTask> { new DbForTask(myDatabase, false) },
                        outputFolder);
                engineToml.Run();

                _bottomUpQValue = new SearchTaskTestCase()
                {
                    Identifier = "BottomUpQValue",
                    OutputDirectory = outputFolder,
                    SpectraFileCount = 2,
                    IsTopDown = false,
                    WriteIndividualResults = false,
                    WritePepXml = false,
                    WritePercolatorFiles = false
                };

                return _bottomUpQValue;
            }
        }

        private static SearchTaskTestCase _bottomUpQValueSingle;
        internal static SearchTaskTestCase BottomUpQValueSingle
        {

            get
            {
                if (!_bottomUpQValueSingle.IsDefault()) return _bottomUpQValueSingle;

                var outputFolder = Path.Combine(ResultDirectory, @"BottomUpQValueSingle");
                string myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                    @"TestData\Task1-SearchTaskconfig.toml");
                SearchTask searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
                string myFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory,
                    @"TestData\TaGe_SA_A549_3_snip_2.mzML");
                string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory,
                    @"TestData\TaGe_SA_A549_3_snip.fasta");
                if (Directory.Exists(outputFolder))
                    Directory.Delete(outputFolder, true);
                var engineToml = new EverythingRunnerEngine(
                    new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) },
                    new List<string> { myFile2 },
                    new List<DbForTask> { new DbForTask(myDatabase, false) },
                    outputFolder);
                engineToml.Run();

                _bottomUpQValueSingle = new SearchTaskTestCase()
                {
                    Identifier = "BottomUpQValueSingle",
                    OutputDirectory = outputFolder,
                    SpectraFileCount = 1,
                    IsTopDown = false,
                    WriteIndividualResults = false,
                    WritePepXml = false,
                    WritePercolatorFiles = false
                };

                return _bottomUpQValueSingle;
            }

        }

        private static SearchTaskTestCase _bottomUpPepQValue;
        internal static SearchTaskTestCase BottomUpPepQValue
        {
            get
            {

                if (!_bottomUpPepQValue.IsDefault()) return _bottomUpPepQValue;

                var outputFolder = Path.Combine(ResultDirectory, @"BottomUpPepQValue");
                if (Directory.Exists(outputFolder))
                    Directory.Delete(outputFolder, true);
                var myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                    @"TestData\Task2-SearchTaskconfig.toml");
                var searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
                string myFile1 = Path.Combine(TestContext.CurrentContext.TestDirectory,
                    @"TestData\TaGe_SA_A549_3_snip.mzML");
                string myFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory,
                    @"TestData\TaGe_SA_A549_3_snip_2.mzML");
                string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory,
                    @"TestData\TaGe_SA_A549_3_snip.fasta");
                // TODO: Uncomment this line and change values for PR 2394
                //searchTaskLoaded.CommonParameters.QValueCutoffForPepCalculation = 0.01;
                var engineToml = new EverythingRunnerEngine(
                    new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) },
                    new List<string> { myFile1, myFile2 }, new List<DbForTask> { new DbForTask(myDatabase, false) },
                    outputFolder);
                engineToml.Run();

                _bottomUpPepQValue = new SearchTaskTestCase()
                {
                    Identifier = "BottomUpPepQValue",
                    OutputDirectory = outputFolder,
                    SpectraFileCount = 2,
                    IsTopDown = false,
                    WriteIndividualResults = false,
                    WritePepXml = false,
                    WritePercolatorFiles = false
                };

                return _bottomUpPepQValue;
            }

        }


        #endregion


    }

    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class PostSearchAnalysisTaskTests
    {
        
        [OneTimeSetUp]
        public static void OneTimeSetUp()
        {
            GlobalVariables.SetUpGlobalVariables();
            if (Directory.Exists(SearchTaskTestCase.ResultDirectory))
                Directory.Delete(SearchTaskTestCase.ResultDirectory, true);
            if (!Directory.Exists(SearchTaskTestCase.ResultDirectory))
                Directory.CreateDirectory(SearchTaskTestCase.ResultDirectory);
        }

        [OneTimeTearDown]
        public static void OneTimeTearDown()
        {
            if (Directory.Exists(SearchTaskTestCase.ResultDirectory))
                Directory.Delete(SearchTaskTestCase.ResultDirectory, true);
        }

        


        /// <summary>
        /// Ensures that there is the proper ratio of summary and individual lines in the result.txt file and that peptides and proteoforms are distinct
        /// </summary>
        /// <param name="allResultTxtLines"></param>
        /// <param name="spectraFileCount"></param>
        /// <param name="isTopDown"></param>
        /// <returns></returns>
        private static string CorrectNumberOfResultLines(string[] allResultTxtLines, int spectraFileCount, bool isTopDown)
        {
            int expectedIndividualFileLines = spectraFileCount == 1 ? 0 : spectraFileCount;
            int expectedSummaryLines = 1;

            var summaryPsmLines = allResultTxtLines.Where(p => p.Contains("All target PSMs")).ToArray();
            var individualPsmLines = allResultTxtLines.Where(p => p.Contains("Target PSMs") 
                                                                  && !p.Contains("All")).ToArray();
            if (summaryPsmLines.Length != expectedSummaryLines)
                return $"Summary Psm Line Mismatch. Was {summaryPsmLines.Length}, Expected {expectedSummaryLines}";
            if (individualPsmLines.Length != expectedIndividualFileLines)
                return $"Individual Psm Line Mismatch. Was {individualPsmLines.Length}, Expected {expectedIndividualFileLines}";

            if (isTopDown)
            {
                var summaryProteoformLines =
                    allResultTxtLines.Where(p => p.Contains("All target proteoforms")).ToArray();
                var individualProteoformLines = allResultTxtLines.Where(p => p.Contains("Target proteoforms")
                                                                             && !p.Contains("All")).ToArray();
                if (summaryProteoformLines.Length != expectedSummaryLines)
                    return $"Summary Proteoform Line Mismatch. Was {summaryProteoformLines.Length}, Expected {expectedSummaryLines}";
                if (individualProteoformLines.Length != expectedIndividualFileLines)
                    return $"Individual Proteoform Line Mismatch. Was {individualProteoformLines.Length}, Expected {expectedIndividualFileLines}";
            }
            else
            {
                var summaryPeptideLines = allResultTxtLines.Where(p => p.Contains("All target peptides")).ToArray();
                var individualPeptideLines = allResultTxtLines.Where(p => p.Contains("Target peptides")
                                                                          && !p.Contains("All")).ToArray();
                if (summaryPeptideLines.Length != expectedSummaryLines)
                    return $"Summary Peptide Line Mismatch. Was {summaryPeptideLines.Length}, Expected {expectedSummaryLines}";
                if (individualPeptideLines.Length != expectedIndividualFileLines)
                    return $"Individual Peptide Line Mismatch. Was {individualPeptideLines.Length}, Expected {expectedIndividualFileLines}";
            }

            var summaryProteinLines = allResultTxtLines.Where(p => p.Contains("All target protein groups")).ToArray();
            var individualProteinLines = allResultTxtLines.Where(p => p.Contains("Target protein groups")
                                                                      && !p.Contains("All")).ToArray();
            if (summaryProteinLines.Length != expectedSummaryLines)
                return $"Summary Protein Line Mismatch. Was {summaryProteinLines.Length}, Expected {expectedSummaryLines}";
            if (individualProteinLines.Length != expectedIndividualFileLines)
                return $"Individual Protein Line Mismatch. Was {individualProteinLines.Length}, Expected {expectedIndividualFileLines}";

            return null;
        }

        [Test]
        public static void AllResultsAndResultsTxtContainsCorrectValues_QValue_BottomUp()
        {
            //First test that AllResults and Results display correct numbers of peptides and psms with q-value filter on
            var outputFolder = SearchTaskTestCase.BottomUpQValue.OutputDirectory;
            string allResultsFile = Path.Combine(outputFolder, "allResults.txt");
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


            string resultsFile = Path.Combine(outputFolder, @"postSearchAnalysisTaskTestOutput\results.txt");
            string[] results = File.ReadAllLines(resultsFile);

            Assert.AreEqual("All target PSMs with q-value <= 0.01: 428", results[5]);
            Assert.AreEqual("All target peptides with q-value <= 0.01: 174", results[6]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target PSMs with q-value <= 0.01: 214", results[9]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target peptides with q-value <= 0.01: 174", results[10]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target protein groups within 1 % FDR: 165", results[11]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target PSMs with q-value <= 0.01: 214", results[13]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target peptides with q-value <= 0.01: 174", results[14]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target protein groups within 1 % FDR: 165", results[15]);
            
            // Search TaGe_SA_A549_3_snip_2 by itself. The results from this should be identical to the file specific results above
            // TaGe_SA_A549_3_snip_2 is searched twice. First with two files being searched simultaneously, then with TaGe_SA_A549_3_snip_2 by itself
            // This allows us to compare the file specific results produced by in the two file search to the output
            // produced by searching the file by itself. The number of PSMs and Peptides observed should be the same
            // for both single-file and multi-file searches.
            // The number of protein groups will be different, because protein inference is performed once, using every peptide
            // identified across all files.
            int TaGe_SA_A549_3_snip_2ExpectedPsms = 214;
            int TaGe_SA_A549_3_snip_2ExpectedPeptides = 174;

            outputFolder = SearchTaskTestCase.BottomUpQValueSingle.OutputDirectory;
            resultsFile = Path.Combine(outputFolder, @"postSearchAnalysisTaskTestOutput\results.txt");
            string[] singleFileResults = File.ReadAllLines(resultsFile);
            Assert.AreEqual("All target PSMs with q-value <= 0.01: " + TaGe_SA_A549_3_snip_2ExpectedPsms, singleFileResults[5]);
            Assert.AreEqual("All target peptides with q-value <= 0.01: " + TaGe_SA_A549_3_snip_2ExpectedPeptides, singleFileResults[6]);
            Assert.AreEqual("All target protein groups with q-value <= 0.01 (1% FDR): 165", singleFileResults[7]);
        }

        [Test]
        public static void AllResultsAndResultsTxtContainsCorrectValues_PepQValue_BottomUp()
        {
            //First test that AllResults and Results display correct numbers of peptides and psms with pep q-value filter on
            var outputFolder = SearchTaskTestCase.BottomUpPepQValue.OutputDirectory;
            var allResultsFile = Path.Combine(outputFolder, "allResults.txt");
            var allResults = File.ReadAllLines(allResultsFile);
            Assert.AreEqual(CorrectNumberOfResultLines(allResults, 2, false), null);
            Assert.AreEqual("All target PSMs with pep q-value <= 0.01: 420", allResults[10]);
            Assert.AreEqual("All target peptides with pep q-value <= 0.01: 172", allResults[11]);
            Assert.AreEqual("All target protein groups with q-value <= 0.01 (1% FDR): 155", allResults[12]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target PSMs with pep q-value <= 0.01: 210", allResults[14]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target peptides with pep q-value <= 0.01: 172", allResults[15]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target protein groups within 1 % FDR: 155", allResults[16]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target PSMs with pep q-value <= 0.01: 210", allResults[18]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target peptides with pep q-value <= 0.01: 172", allResults[19]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target protein groups within 1 % FDR: 155", allResults[20]);

            var resultsFile = Path.Combine(outputFolder, @"postSearchAnalysisTaskTestOutput\results.txt");
            var results = File.ReadAllLines(resultsFile);
            Assert.AreEqual(CorrectNumberOfResultLines(results, 2, false), null);
            Assert.AreEqual("All target PSMs with pep q-value <= 0.01: 420", results[5]);
            Assert.AreEqual("All target peptides with pep q-value <= 0.01: 172", results[6]);
            Assert.AreEqual("All target protein groups with q-value <= 0.01 (1% FDR): 155", results[7]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target PSMs with pep q-value <= 0.01: 210", results[9]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target peptides with pep q-value <= 0.01: 172", results[10]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target protein groups within 1 % FDR: 155", results[11]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target PSMs with pep q-value <= 0.01: 210", results[13]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target peptides with pep q-value <= 0.01: 172", results[14]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target protein groups within 1 % FDR: 155", results[15]);
        }


        [Test]
        [TestCaseSource(typeof(SearchTaskTestCase), nameof(SearchTaskTestCase.AllTestCaseIdentifiers))]
        public static void AllResultAndResultTxtContainsCorrectNumberOfLines(string testCaseIdentifier)
        {
            var testCase = SearchTaskTestCase.TestCaseLocator(testCaseIdentifier);
        }
    }
}