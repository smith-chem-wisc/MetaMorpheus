using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework; using Assert = NUnit.Framework.Legacy.ClassicAssert;

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

            EverythingRunnerEngineTestCase.TryGetTestCase(EverythingRunnerEngineTestCases.BottomUpQValueSingle, out var testCaseSingle);
            outputFolder = testCaseSingle.OutputDirectory;

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
            EverythingRunnerEngineTestCase.TryGetTestCase(EverythingRunnerEngineTestCases.BottomUpPepQValue, out var testCase);
            string outputFolder = testCase.OutputDirectory;
            var allResultsFile = Path.Combine(outputFolder, "allResults.txt");
            var allResults = File.ReadAllLines(allResultsFile);
            Assert.AreEqual("All target PSMs with pep q-value <= 0.01: 382", allResults[10]);
            Assert.AreEqual("All target peptides with pep q-value <= 0.01: 153", allResults[11]);
            Assert.AreEqual("All target protein groups with q-value <= 0.01 (1% FDR): 140", allResults[12]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target PSMs with pep q-value <= 0.01: 190", allResults[14]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target peptides with pep q-value <= 0.01: 153", allResults[15]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target protein groups within 1 % FDR: 140", allResults[16]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target PSMs with pep q-value <= 0.01: 190", allResults[18]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target peptides with pep q-value <= 0.01: 153", allResults[19]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target protein groups within 1 % FDR: 140", allResults[20]);

            var resultsFile = Path.Combine(outputFolder, @"postSearchAnalysisTaskTestOutput\results.txt");
            var results = File.ReadAllLines(resultsFile);
            Assert.AreEqual("All target PSMs with pep q-value <= 0.01: 382", results[5]);
            Assert.AreEqual("All target peptides with pep q-value <= 0.01: 153", results[6]);
            Assert.AreEqual("All target protein groups with q-value <= 0.01 (1% FDR): 140", results[7]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target PSMs with pep q-value <= 0.01: 190", results[9]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target peptides with pep q-value <= 0.01: 153", results[10]);
            Assert.AreEqual("TaGe_SA_A549_3_snip - Target protein groups within 1 % FDR: 140", results[11]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target PSMs with pep q-value <= 0.01: 190", results[13]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target peptides with pep q-value <= 0.01: 153", results[14]);
            Assert.AreEqual("TaGe_SA_A549_3_snip_2 - Target protein groups within 1 % FDR: 140", results[15]);
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
            Assert.AreEqual(expectedSummaryLines, summaryPsmLines.Length);
            Assert.AreEqual(expectedIndividualFileLines, individualPsmLines.Length);

            if (testCase.IsTopDown)
            {
                var summaryProteoformLines =
                    allResultTxtLines.Where(p => p.Contains("All target proteoforms")).ToArray();
                var individualProteoformLines = allResultTxtLines.Where(p => p.Contains("Target proteoforms")
                                                                             && !p.Contains("All")).ToArray();
                Assert.AreEqual(expectedSummaryLines, summaryProteoformLines.Length);
                Assert.AreEqual(expectedIndividualFileLines, individualProteoformLines.Length);
            }
            else
            {
                var summaryPeptideLines = allResultTxtLines.Where(p => p.Contains("All target peptides")).ToArray();
                var individualPeptideLines = allResultTxtLines.Where(p => p.Contains("Target peptides")
                                                                          && !p.Contains("All")).ToArray();
                Assert.AreEqual(expectedSummaryLines, summaryPeptideLines.Length);
                Assert.AreEqual(expectedIndividualFileLines, individualPeptideLines.Length);
            }

            var summaryProteinLines = allResultTxtLines.Where(p => p.Contains("All target protein groups")).ToArray();
            var individualProteinLines = allResultTxtLines.Where(p => p.Contains("Target protein groups")
                                                                      && !p.Contains("All")).ToArray();
            Assert.AreEqual(expectedSummaryLines, summaryProteinLines.Length);
            Assert.AreEqual(expectedIndividualFileLines, individualProteinLines.Length);
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

            Assert.AreEqual(expectedResultFileCount, psmFiles.Length);
            Assert.AreEqual(expectedResultFileCount, proteinGroupFiles.Length);
            if (testCase.IsTopDown)
            {
                Assert.AreEqual(expectedResultFileCount, proteoformFiles.Length);
                Assert.AreEqual(0, peptideFiles.Length);
            }
            else
            {
                Assert.AreEqual(expectedResultFileCount, peptideFiles.Length);
                Assert.AreEqual(0, proteoformFiles.Length);
            }

            if (testCase.WritePepXml)
            {
                Assert.AreEqual(spectraFileCount, pepXmlFiles.Length);
            }
            else
            {
                Assert.AreEqual(0, pepXmlFiles.Length);
            }

            if (testCase.WriteIndividualResults)
            {
                Assert.AreEqual(expectedResultFileCount, percolatorFiles.Length);
            }
            else
            {
                Assert.AreEqual(1, percolatorFiles.Length);
            }

            if (testCase.WriteMzId)
            {
                Assert.AreEqual(spectraFileCount, mzidFiles.Length);
            }
            else
            {
                Assert.AreEqual(0, mzidFiles.Length);
            }
        }
    }
}