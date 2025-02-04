using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Uses test cases found in EverythingRunnerEngineTestCase.cs
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class PostSearchAnalysisTaskTests
    {
        public static Array GetTestCases() => Enum.GetValues(typeof(EverythingRunnerEngineTestCases));

        [Test]
        public static void AllResultsAndResultsTxtContainsCorrectValues_QValue_BottomUp()
        {
            //First test that AllResults and Results display correct numbers of peptides and psms with q-value filter on
            EverythingRunnerEngineTestCase.TryGetTestCase(EverythingRunnerEngineTestCases.BottomUpQValue, out var testCase);
            string outputFolder = testCase.OutputDirectory;
            string allResultsFile = Path.Combine(outputFolder, "allResults.txt");
            string[] allResults = File.ReadAllLines(allResultsFile);

            // The new PEP calculation method improves things, so all these numbers are increasing as of (7/17/24)
            // There is a discrepancy between the number of All target peptides and individual file target peptides, 
            // presumably due to the way protein inference is performed.
            Assert.That(allResults[10], Is.EqualTo("All target PSMs with q-value <= 0.01: 428"));
            Assert.That(allResults[11], Is.EqualTo("All target peptides with q-value <= 0.01: 174"));
            Assert.That(allResults[12], Is.EqualTo("All target protein groups with q-value <= 0.01 (1% FDR): 165"));
            Assert.That(allResults[14], Is.EqualTo("TaGe_SA_A549_3_snip - Target PSMs with q-value <= 0.01: 214"));
            Assert.That(allResults[15], Is.EqualTo("TaGe_SA_A549_3_snip - Target peptides with q-value <= 0.01: 174"));
            Assert.That(allResults[16], Is.EqualTo("TaGe_SA_A549_3_snip - Target protein groups within 1 % FDR: 165"));
            Assert.That(allResults[18], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target PSMs with q-value <= 0.01: 214"));
            Assert.That(allResults[19], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target peptides with q-value <= 0.01: 174"));
            Assert.That(allResults[20], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target protein groups within 1 % FDR: 165"));


            string resultsFile = Path.Combine(outputFolder, @"postSearchAnalysisTaskTestOutput\results.txt");
            string[] results = File.ReadAllLines(resultsFile);

            Assert.That(results[5], Is.EqualTo("All target PSMs with q-value <= 0.01: 428"));
            Assert.That(results[6], Is.EqualTo("All target peptides with q-value <= 0.01: 174"));
            Assert.That(results[9], Is.EqualTo("TaGe_SA_A549_3_snip - Target PSMs with q-value <= 0.01: 214"));
            Assert.That(results[10], Is.EqualTo("TaGe_SA_A549_3_snip - Target peptides with q-value <= 0.01: 174"));
            Assert.That(results[11], Is.EqualTo("TaGe_SA_A549_3_snip - Target protein groups within 1 % FDR: 165"));
            Assert.That(results[13], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target PSMs with q-value <= 0.01: 214"));
            Assert.That(results[14], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target peptides with q-value <= 0.01: 174"));
            Assert.That(results[15], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target protein groups within 1 % FDR: 165"));

            // Search TaGe_SA_A549_3_snip_2 by itself. The results from this should be identical to the file specific results above
            // TaGe_SA_A549_3_snip_2 is searched twice. First with two files being searched simultaneously, then with TaGe_SA_A549_3_snip_2 by itself
            // This allows us to compare the file specific results produced by in the two file search to the output
            // produced by searching the file by itself. The number of PSMs and Peptides observed should be the same
            // for both single-file and multi-file searches.
            // The number of protein groups will be different, because protein inference is performed once, using every peptide
            // identified across all files.
            int TaGe_SA_A549_3_snip_2ExpectedPsms = 214;
            int TaGe_SA_A549_3_snip_2ExpectedPeptides = 174;

            EverythingRunnerEngineTestCase.TryGetTestCase(EverythingRunnerEngineTestCases.BottomUpQValueSingle, out var testCaseSingle);
            outputFolder = testCaseSingle.OutputDirectory;

            resultsFile = Path.Combine(outputFolder, @"postSearchAnalysisTaskTestOutput\results.txt");
            string[] singleFileResults = File.ReadAllLines(resultsFile);
            Assert.That(singleFileResults[5], Is.EqualTo("All target PSMs with q-value <= 0.01: " + TaGe_SA_A549_3_snip_2ExpectedPsms));
            Assert.That(singleFileResults[6], Is.EqualTo("All target peptides with q-value <= 0.01: " + TaGe_SA_A549_3_snip_2ExpectedPeptides));
            Assert.That(singleFileResults[7], Is.EqualTo("All target protein groups with q-value <= 0.01 (1% FDR): 165"));
        }

        [Test]
        public static void LocalTest()
        {
            var modelPath = @"D:\MetaMorpheusVignette\Search_0105GPTMD_NotchSorting\Task1-SearchTask\model.zip";
            //mLContext.Model.Load(Path.Combine(outputFolder, "model.zip"), out DataViewSchema savedModelSchema);
        }

        [Test]
        public static void AllResultsAndResultsTxtContainsCorrectValues_PepQValue_BottomUp()
        {
            //First test that AllResults and Results display correct numbers of peptides and psms with pep q-value filter on
            EverythingRunnerEngineTestCase.TryGetTestCase(EverythingRunnerEngineTestCases.BottomUpPepQValue, out var testCase);
            string outputFolder = testCase.OutputDirectory;
            var allResultsFile = Path.Combine(outputFolder, "allResults.txt");
            var allResults = File.ReadAllLines(allResultsFile);
            Assert.That(allResults[10], Is.EqualTo("All target PSMs with pep q-value <= 0.01: 382"));
            Assert.That(allResults[11], Is.EqualTo("All target peptides with pep q-value <= 0.01: 153"));
            Assert.That(allResults[12], Is.EqualTo("All target protein groups with q-value <= 0.01 (1% FDR): 140"));
            Assert.That(allResults[14], Is.EqualTo("TaGe_SA_A549_3_snip - Target PSMs with pep q-value <= 0.01: 190"));
            Assert.That(allResults[15], Is.EqualTo("TaGe_SA_A549_3_snip - Target peptides with pep q-value <= 0.01: 153"));
            Assert.That(allResults[16], Is.EqualTo("TaGe_SA_A549_3_snip - Target protein groups within 1 % FDR: 140"));
            Assert.That(allResults[18], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target PSMs with pep q-value <= 0.01: 190"));
            Assert.That(allResults[19], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target peptides with pep q-value <= 0.01: 153"));
            Assert.That(allResults[20], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target protein groups within 1 % FDR: 140"));

            var resultsFile = Path.Combine(outputFolder, @"postSearchAnalysisTaskTestOutput\results.txt");
            var results = File.ReadAllLines(resultsFile);
            Assert.That(results[5], Is.EqualTo("All target PSMs with pep q-value <= 0.01: 382"));
            Assert.That(results[6], Is.EqualTo("All target peptides with pep q-value <= 0.01: 153"));
            Assert.That(results[7], Is.EqualTo("All target protein groups with q-value <= 0.01 (1% FDR): 140"));
            Assert.That(results[9], Is.EqualTo("TaGe_SA_A549_3_snip - Target PSMs with pep q-value <= 0.01: 190"));
            Assert.That(results[10], Is.EqualTo("TaGe_SA_A549_3_snip - Target peptides with pep q-value <= 0.01: 153"));
            Assert.That(results[11], Is.EqualTo("TaGe_SA_A549_3_snip - Target protein groups within 1 % FDR: 140"));
            Assert.That(results[13], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target PSMs with pep q-value <= 0.01: 190"));
            Assert.That(results[14], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target peptides with pep q-value <= 0.01: 153"));
            Assert.That(results[15], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target protein groups within 1 % FDR: 140"));
        }

        [Test]
        public static void PepAnalysisEngineHasReproducibleOutput()
        {
            EverythingRunnerEngineTestCase.TryGetTestCase(EverythingRunnerEngineTestCases.BottomUpPepQValue, out var testCase);
            string outputFolderPepQ = testCase.OutputDirectory;
            var resultsFilePepQ = Path.Combine(outputFolderPepQ, @"postSearchAnalysisTaskTestOutput\results.txt");
            string[] resultsPepQ = File.ReadAllLines(resultsFilePepQ);

            //First test that AllResults and Results display correct numbers of peptides and psms with q-value filter on
            EverythingRunnerEngineTestCase.TryGetTestCase(EverythingRunnerEngineTestCases.BottomUpQValue, out testCase);
            string outputFolderQ = testCase.OutputDirectory;
            string resultsFileQ = Path.Combine(outputFolderQ, @"postSearchAnalysisTaskTestOutput\results.txt");
            string[] resultsQ = File.ReadAllLines(resultsFileQ);

            // The Q-Value test case and the PEP-Q value test case represent the same data ran twice
            // These assert staments compare the outputs of the PepAnalysisEngine for each test case
            // They should be identical!!! If they aren't, then PEP is not reproducible for some reason
            Assert.That(resultsPepQ[37], Is.EqualTo(resultsQ[37]));
            Assert.That(resultsPepQ[38], Is.EqualTo(resultsQ[38]));
            Assert.That(resultsPepQ[39], Is.EqualTo(resultsQ[39]));
            Assert.That(resultsPepQ[40], Is.EqualTo(resultsQ[40]));
            Assert.That(resultsPepQ[41], Is.EqualTo(resultsQ[41]));
            Assert.That(resultsPepQ[42], Is.EqualTo(resultsQ[42]));
            Assert.That(resultsPepQ[43], Is.EqualTo(resultsQ[43]));
            Assert.That(resultsPepQ[44], Is.EqualTo(resultsQ[44]));
            Assert.That(resultsPepQ[45], Is.EqualTo(resultsQ[45]));
            Assert.That(resultsPepQ[46], Is.EqualTo(resultsQ[46]));
            Assert.That(resultsPepQ[47], Is.EqualTo(resultsQ[47]));
            Assert.That(resultsPepQ[48], Is.EqualTo(resultsQ[48]));
            Assert.That(resultsPepQ[49], Is.EqualTo(resultsQ[49]));
            Assert.That(resultsPepQ[50], Is.EqualTo(resultsQ[50]));
        }

        /// <summary>
        /// Ensures that there is the proper ratio of summary and individual lines in the result.txt file and that peptides and proteoforms are distinct
        /// </summary>
        [Test]
        [TestCaseSource(nameof(GetTestCases))]
        public static void AllResultTxtContainsCorrectNumberOfResultLines(EverythingRunnerEngineTestCases testCaseIdentifier)
        {
            var testCase = EverythingRunnerEngineTestCase.GetTestCase(testCaseIdentifier);

            int expectedIndividualFileLines = testCase.DataFileList.Count == 1 || !testCase.WriteIndividualResults
                ? 0 : testCase.DataFileList.Count;
            int expectedSummaryLines = 1;
            var allResultTxtLines = File.ReadAllLines(Path.Combine(testCase.OutputDirectory, @"allResults.txt"));

            var summaryPsmLines = allResultTxtLines.Where(p => p.Contains("All target PSMs")).ToArray();
            var individualPsmLines = allResultTxtLines.Where(p => p.Contains("Target PSMs")
                                                                  && !p.Contains("All")).ToArray();
            Assert.That(summaryPsmLines.Length, Is.EqualTo(expectedSummaryLines));
            Assert.That(individualPsmLines.Length, Is.EqualTo(expectedIndividualFileLines));

            if (testCase.IsTopDown)
            {
                var summaryProteoformLines =
                    allResultTxtLines.Where(p => p.Contains("All target proteoforms")).ToArray();
                var individualProteoformLines = allResultTxtLines.Where(p => p.Contains("Target proteoforms")
                                                                             && !p.Contains("All")).ToArray();
                Assert.That(summaryProteoformLines.Length, Is.EqualTo(expectedSummaryLines));
                Assert.That(individualProteoformLines.Length, Is.EqualTo(expectedIndividualFileLines));
            }
            else
            {
                var summaryPeptideLines = allResultTxtLines.Where(p => p.Contains("All target peptides")).ToArray();
                var individualPeptideLines = allResultTxtLines.Where(p => p.Contains("Target peptides")
                                                                          && !p.Contains("All")).ToArray();
                Assert.That(summaryPeptideLines.Length, Is.EqualTo(expectedSummaryLines));
                Assert.That(individualPeptideLines.Length, Is.EqualTo(expectedIndividualFileLines));
            }

            var summaryProteinLines = allResultTxtLines.Where(p => p.Contains("All target protein groups")).ToArray();
            var individualProteinLines = allResultTxtLines.Where(p => p.Contains("Target protein groups")
                                                                      && !p.Contains("All")).ToArray();
            Assert.That(summaryProteinLines.Length, Is.EqualTo(expectedSummaryLines));
            Assert.That(individualProteinLines.Length, Is.EqualTo(expectedIndividualFileLines));
        }

        /// <summary>
        /// Ensures that the files written out with each search are correct according to the search parameters and data type
        /// </summary>
        [Test]
        [TestCaseSource(nameof(GetTestCases))]
        public static void CorrectFilesAreWrittenWithCorrectName(EverythingRunnerEngineTestCases testCaseIdentifier)
        {
            var testCase = EverythingRunnerEngineTestCase.GetTestCase(testCaseIdentifier);
            var psmFiles = Directory.GetFiles(testCase.OutputDirectory, "*PSMs.psmtsv", SearchOption.AllDirectories);
            var pepXmlFiles = Directory.GetFiles(testCase.OutputDirectory, "*.pep.xml", SearchOption.AllDirectories);
            var percolatorFiles = Directory.GetFiles(testCase.OutputDirectory, "*Percolator.tab", SearchOption.AllDirectories);
            var proteinGroupFiles = Directory.GetFiles(testCase.OutputDirectory, "*ProteinGroups.tsv", SearchOption.AllDirectories);
            var peptideFiles = Directory.GetFiles(testCase.OutputDirectory, "*Peptides.psmtsv", SearchOption.AllDirectories);
            var proteoformFiles = Directory.GetFiles(testCase.OutputDirectory, "*Proteoforms.psmtsv", SearchOption.AllDirectories);
            var mzidFiles = Directory.GetFiles(testCase.OutputDirectory, "*.mzid", SearchOption.AllDirectories);

            int spectraFileCount = testCase.DataFileList.Count;
            var expectedResultFileCount = testCase.WriteIndividualResults && testCase.DataFileList.Count > 1
                ? testCase.DataFileList.Count + 1 : 1;

            Assert.That(psmFiles.Length, Is.EqualTo(expectedResultFileCount));
            Assert.That(proteinGroupFiles.Length, Is.EqualTo(expectedResultFileCount));
            if (testCase.IsTopDown)
            {
                Assert.That(proteoformFiles.Length, Is.EqualTo(expectedResultFileCount));
                Assert.That(peptideFiles.Length, Is.EqualTo(0));
            }
            else
            {
                Assert.That(peptideFiles.Length, Is.EqualTo(expectedResultFileCount));
                Assert.That(proteoformFiles.Length, Is.EqualTo(0));
            }

            if (testCase.WritePepXml)
            {
                Assert.That(pepXmlFiles.Length, Is.EqualTo(spectraFileCount));
            }
            else
            {
                Assert.That(pepXmlFiles.Length, Is.EqualTo(0));
            }

            if (testCase.WriteIndividualResults)
            {
                Assert.That(percolatorFiles.Length, Is.EqualTo(expectedResultFileCount));
            }
            else
            {
                Assert.That(percolatorFiles.Length, Is.EqualTo(1));
            }

            if (testCase.WriteMzId)
            {
                Assert.That(mzidFiles.Length, Is.EqualTo(spectraFileCount));
            }
            else
            {
                Assert.That(mzidFiles.Length, Is.EqualTo(0));
            }
        }
    }
}