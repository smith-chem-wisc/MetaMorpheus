using Easy.Common.Extensions;
using Chemistry;
using EngineLayer;
using EngineLayer.SpectrumMatch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Reflection;
using EngineLayer.DatabaseLoading;
using TaskLayer;
using UsefulProteomicsDatabases;

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

            Assert.That(allResults[10], Is.EqualTo("All target PSMs with q-value <= 0.01: 428"));
            Assert.That(allResults[11], Is.EqualTo("All target peptides with q-value <= 0.01: 173"));
            Assert.That(allResults[12], Is.EqualTo("All target protein groups with q-value <= 0.01 (1% FDR): 165"));
            Assert.That(allResults[13], Is.EqualTo("All Precursors: 1070"));
            Assert.That(allResults[14], Is.EqualTo("All MS2 Scans: 294"));
            Assert.That(allResults[16], Is.EqualTo("TaGe_SA_A549_3_snip - Target PSMs with q-value <= 0.01: 214"));
            Assert.That(allResults[17], Is.EqualTo("TaGe_SA_A549_3_snip - Target peptides with q-value <= 0.01: 174"));
            Assert.That(allResults[18], Is.EqualTo("TaGe_SA_A549_3_snip - Target protein groups with q-value <= 0.01: 165"));
            Assert.That(allResults[19], Is.EqualTo("TaGe_SA_A549_3_snip - Precursors: 535"));
            Assert.That(allResults[20], Is.EqualTo("TaGe_SA_A549_3_snip - MS2 Scans: 147"));
            Assert.That(allResults[22], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target PSMs with q-value <= 0.01: 214"));
            Assert.That(allResults[23], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target peptides with q-value <= 0.01: 174"));
            Assert.That(allResults[24], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target protein groups with q-value <= 0.01: 165"));
            Assert.That(allResults[25], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Precursors: 535"));
            Assert.That(allResults[26], Is.EqualTo("TaGe_SA_A549_3_snip_2 - MS2 Scans: 147"));

            string resultsFile = Path.Combine(outputFolder, @"postSearchAnalysisTaskTestOutput\results.txt");
            string[] results = File.ReadAllLines(resultsFile);

            Assert.That(results[5], Is.EqualTo("All target PSMs with q-value <= 0.01: 428"));
            Assert.That(results[6], Is.EqualTo("All target peptides with q-value <= 0.01: 173"));
            Assert.That(results[7], Is.EqualTo("All target protein groups with q-value <= 0.01 (1% FDR): 165"));
            Assert.That(results[8], Is.EqualTo("All Precursors: 1070"));
            Assert.That(results[9], Is.EqualTo("All MS2 Scans: 294"));

            Assert.That(results[11], Is.EqualTo("TaGe_SA_A549_3_snip - Target PSMs with q-value <= 0.01: 214"));
            Assert.That(results[12], Is.EqualTo("TaGe_SA_A549_3_snip - Target peptides with q-value <= 0.01: 174"));
            Assert.That(results[13], Is.EqualTo("TaGe_SA_A549_3_snip - Target protein groups with q-value <= 0.01: 165"));
            Assert.That(results[14], Is.EqualTo("TaGe_SA_A549_3_snip - Precursors: 535"));
            Assert.That(results[15], Is.EqualTo("TaGe_SA_A549_3_snip - MS2 Scans: 147"));

            Assert.That(results[17], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target PSMs with q-value <= 0.01: 214"));
            Assert.That(results[18], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target peptides with q-value <= 0.01: 174"));
            Assert.That(results[19], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target protein groups with q-value <= 0.01: 165"));
            Assert.That(results[20], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Precursors: 535"));
            Assert.That(results[21], Is.EqualTo("TaGe_SA_A549_3_snip_2 - MS2 Scans: 147"));


            // Search TaGe_SA_A549_3_snip_2 by itself. The results from this should be identical to the file specific results above
            // TaGe_SA_A549_3_snip_2 is searched twice. First with two files being searched simultaneously, then with TaGe_SA_A549_3_snip_2 by itself
            // This allows us to compare the file specific results produced by in the two file search to the output
            // produced by searching the file by itself. The number of PSMs and Peptides observed should be the same
            // for both single-file and multi-file searches.
            // The number of protein groups will be different, because protein inference is performed once, using every peptide
            // identified across all files.
            int TaGe_SA_A549_3_snip_2ExpectedPsms = 214;
            int TaGe_SA_A549_3_snip_2ExpectedPeptides = 173;

            EverythingRunnerEngineTestCase.TryGetTestCase(EverythingRunnerEngineTestCases.BottomUpQValueSingle, out var testCaseSingle);
            outputFolder = testCaseSingle.OutputDirectory;

            resultsFile = Path.Combine(outputFolder, @"postSearchAnalysisTaskTestOutput\results.txt");
            string[] singleFileResults = File.ReadAllLines(resultsFile);
            Assert.That(singleFileResults[5], Is.EqualTo("All target PSMs with q-value <= 0.01: " + TaGe_SA_A549_3_snip_2ExpectedPsms));
            Assert.That(singleFileResults[6], Is.EqualTo("All target peptides with q-value <= 0.01: " + TaGe_SA_A549_3_snip_2ExpectedPeptides));
            Assert.That(singleFileResults[7], Is.EqualTo("All target protein groups with q-value <= 0.01 (1% FDR): 165"));
        }

        [Test]
        public static void AllResultsAndResultsTxtContainsCorrectValues_PepQValue_BottomUp()
        {
            //First test that AllResults and Results display correct numbers of peptides and psms with pep q-value filter on
            // Note: PEP Q-value filtering typically yields fewer PSMs, peptides, and protein groups than Q-value filtering
            // (e.g., 144 protein groups here vs. 165 with Q-value filtering) because PEP is a more stringent per-PSM
            // confidence metric that considers additional features beyond just target-decoy competition.
            EverythingRunnerEngineTestCase.TryGetTestCase(EverythingRunnerEngineTestCases.BottomUpPepQValue, out var testCase);
            string outputFolder = testCase.OutputDirectory;
            var allResultsFile = Path.Combine(outputFolder, "allResults.txt");
            var allResults = File.ReadAllLines(allResultsFile);

            Assert.That(allResults[10], Is.EqualTo("All target PSMs with pep q-value <= 0.01: 374"));
            Assert.That(allResults[11], Is.EqualTo("All target peptides with pep q-value <= 0.01: 153"));
            Assert.That(allResults[12], Is.EqualTo("All target protein groups with pep q-value <= 0.01: 144"));
            Assert.That(allResults[13], Is.EqualTo("All Precursors: 1070"));
            Assert.That(allResults[14], Is.EqualTo("All MS2 Scans: 294"));

            Assert.That(allResults[16], Is.EqualTo("TaGe_SA_A549_3_snip - Target PSMs with pep q-value <= 0.01: 187"));
            Assert.That(allResults[17], Is.EqualTo("TaGe_SA_A549_3_snip - Target peptides with pep q-value <= 0.01: 151"));
            Assert.That(allResults[18], Is.EqualTo("TaGe_SA_A549_3_snip - Target protein groups with pep q-value <= 0.01: 144"));
            Assert.That(allResults[19], Is.EqualTo("TaGe_SA_A549_3_snip - Precursors: 535"));
            Assert.That(allResults[20], Is.EqualTo("TaGe_SA_A549_3_snip - MS2 Scans: 147"));

            Assert.That(allResults[22], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target PSMs with pep q-value <= 0.01: 187"));
            Assert.That(allResults[23], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target peptides with pep q-value <= 0.01: 151"));
            Assert.That(allResults[24], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target protein groups with pep q-value <= 0.01: 144"));
            Assert.That(allResults[25], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Precursors: 535"));
            Assert.That(allResults[26], Is.EqualTo("TaGe_SA_A549_3_snip_2 - MS2 Scans: 147"));


            var resultsFile = Path.Combine(outputFolder, @"postSearchAnalysisTaskTestOutput\results.txt");
            var results = File.ReadAllLines(resultsFile);

            Assert.That(results[5], Is.EqualTo("All target PSMs with pep q-value <= 0.01: 374"));
            Assert.That(results[6], Is.EqualTo("All target peptides with pep q-value <= 0.01: 153"));
            Assert.That(results[7], Is.EqualTo("All target protein groups with pep q-value <= 0.01: 144"));
            Assert.That(results[8], Is.EqualTo("All Precursors: 1070"));
            Assert.That(results[9], Is.EqualTo("All MS2 Scans: 294"));

            Assert.That(results[11], Is.EqualTo("TaGe_SA_A549_3_snip - Target PSMs with pep q-value <= 0.01: 187"));
            Assert.That(results[12], Is.EqualTo("TaGe_SA_A549_3_snip - Target peptides with pep q-value <= 0.01: 151"));
            Assert.That(results[13], Is.EqualTo("TaGe_SA_A549_3_snip - Target protein groups with pep q-value <= 0.01: 144"));
            Assert.That(results[14], Is.EqualTo("TaGe_SA_A549_3_snip - Precursors: 535"));
            Assert.That(results[15], Is.EqualTo("TaGe_SA_A549_3_snip - MS2 Scans: 147"));

            Assert.That(results[17], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target PSMs with pep q-value <= 0.01: 187"));
            Assert.That(results[18], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target peptides with pep q-value <= 0.01: 151"));
            Assert.That(results[19], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Target protein groups with pep q-value <= 0.01: 144"));
            Assert.That(results[20], Is.EqualTo("TaGe_SA_A549_3_snip_2 - Precursors: 535"));
            Assert.That(results[21], Is.EqualTo("TaGe_SA_A549_3_snip_2 - MS2 Scans: 147"));
        }

        [Test]
        public static void AllResultsAndResultsTxtWorksForTimsTof()
        {
            //First test that AllResults and Results display correct numbers of peptides and psms with q-value filter on
            EverythingRunnerEngineTestCase.TryGetTestCase(EverythingRunnerEngineTestCases.timsTOFRawFileMixed, out var testCase);
            string outputFolder = testCase.OutputDirectory;
            string allResultsFile = Path.Combine(outputFolder, "allResults.txt");
            string[] allResults = File.ReadAllLines(allResultsFile);

            // The database used for the timsTof file is tiny, only 14 proteins. As a results, we don't get any confident results at the peptide level
            // This is fine for the purposes of this test
            Assert.That(allResults[10], Is.EqualTo("All target PSMs with q-value <= 0.01: 143"));
            Assert.That(allResults[16], Is.EqualTo("snippet - Target PSMs with q-value <= 0.01: 120"));
            Assert.That(allResults[22], Is.EqualTo("TaGe_SA_A549_3_snip - Target PSMs with q-value <= 0.01: 187"));

            string resultsFile = Path.Combine(outputFolder, @"postSearchAnalysisTaskTestOutput\results.txt");
            string[] results = File.ReadAllLines(resultsFile);

            Assert.That(results[5], Is.EqualTo("All target PSMs with q-value <= 0.01: 143"));
            Assert.That(results[11], Is.EqualTo("snippet - Target PSMs with q-value <= 0.01: 120"));
            Assert.That(results[17], Is.EqualTo("TaGe_SA_A549_3_snip - Target PSMs with q-value <= 0.01: 187"));

            string allPsmsFile = Path.Combine(outputFolder, "postSearchAnalysisTaskTestOutput", "AllPSMs.psmtsv");
            string[] allPsms = File.ReadAllLines(allPsmsFile);
            Assert.That(allPsms[0].Split('\t').Any(header => header.Equals("1/K0")));
            int indexOfFirstRawPSM = allPsms.IndexOf(allPsms.First(l => l.StartsWith("TaGe_SA_A549_3")));
            int indexOfFirstTOFPSM = allPsms.IndexOf(allPsms.First(l => l.StartsWith("snippet")));
            int indexOf1overK0Column = allPsms[0].Split('\t').ToList().IndexOf("1/K0");
            Assert.That(allPsms[indexOfFirstRawPSM].Split('\t')[indexOf1overK0Column], Is.EqualTo("N/A"));
            Assert.That(double.Parse(allPsms[indexOfFirstTOFPSM].Split('\t')[indexOf1overK0Column]), Is.EqualTo(1).Within(2)); // Don't actually care about the value, just want to ensure it's a double

            string rawPsmsFile = Path.Combine(outputFolder, "postSearchAnalysisTaskTestOutput", "Individual File Results", "TaGe_SA_A549_3_snip_PSMs.psmtsv");
            string[] rawPsms = File.ReadAllLines(rawPsmsFile);
            Assert.That(!rawPsms[0].Split('\t').Any(header => header.Equals("1/K0"))); // Assert that 1/k0 column doesn't get written unless it's needed
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
            int startIndexQ = resultsQ.IndexOf("Engine type: FdrAnalysisEngine");
            int startIndexPep = resultsPepQ.IndexOf("Engine type: FdrAnalysisEngine");

            Assert.That(startIndexQ, Is.EqualTo(startIndexPep));

            int starsFound = 0;
            for (int i = startIndexQ; i < int.MaxValue; i++)
            {
                var qLine = resultsQ[i];
                var pepLine = resultsPepQ[i];

                if (qLine.StartsWith("Time to run"))
                    continue;

                if (qLine.StartsWith("*****"))
                {
                    starsFound++;
                    if (starsFound == 2)
                    {
                        break; // We've reached the end of the PepAnalysisEngine output
                    }
                }

                Assert.That(pepLine, Is.EqualTo(qLine), "Outputs of PepAnalysisEngine differ between Q-Value and PEP Q-Value test cases");
            }
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

        [Test]
        public static void WriteDigestionCountsByProtein_WritesCorrectFile()
        {
            // Arrange
            var task = new PostSearchAnalysisTask();
            var outputDirectory = Path.Combine(TestContext.CurrentContext.WorkDirectory, "DigestionCountTest");
            if (Directory.Exists(outputDirectory))
                Directory.Delete(outputDirectory, true);
            Directory.CreateDirectory(outputDirectory);
            var parameters = new PostSearchAnalysisParameters
            {
                SearchParameters = new(),
                OutputFolder = outputDirectory,
                SearchTaskId = "TestTask",
            };

            task.GetType().GetProperty("Parameters").SetValue(task, parameters);
            var digestionCountDictionary = new Dictionary<(string Accession, string BaseSeqeunce), int>
            {
                { ("Protein1", "SEQUENCE1"), 5 },
                { ("Protein2", "SEQUENCE2"), 10 }
            };
            task.GetType().GetProperty("DigestionCountDictionary", BindingFlags.NonPublic | BindingFlags.Instance).SetValue(task, digestionCountDictionary);

            // Act
            var method = task.GetType().GetMethod("WriteDigestionCountByProtein", BindingFlags.NonPublic | BindingFlags.Instance);
            method!.Invoke(task, null);

            // Assert
            var expectedFilePath = Path.Combine(parameters.OutputFolder, "DigestionCountsByProteins.tsv");
            Assert.That(File.Exists(expectedFilePath), Is.True);

            var lines = File.ReadAllLines(expectedFilePath);
            Assert.That(lines.Length, Is.EqualTo(3));
            Assert.That(lines[0], Is.EqualTo("Protein Accession\tPrimary Sequence\tDigestion Products"));
            Assert.That(lines[1], Is.EqualTo("Protein1\tSEQUENCE1\t5"));
            Assert.That(lines[2], Is.EqualTo("Protein2\tSEQUENCE2\t10"));

            // Cleanup
            Directory.Delete(parameters.OutputFolder, true);
        }

        [Test]
        public static void WriteDigestionCountsHistogram_WritesCorrectFile()
        {
            // Arrange
            var task = new PostSearchAnalysisTask() { CommonParameters = new() };
            var outputDirectory = Path.Combine(TestContext.CurrentContext.WorkDirectory, "DigestionHistogramTest");
            if (Directory.Exists(outputDirectory))
                Directory.Delete(outputDirectory, true);
            Directory.CreateDirectory(outputDirectory);
            var parameters = new PostSearchAnalysisParameters
            {
                SearchParameters = new(),
                OutputFolder = outputDirectory,
                SearchTaskId = "TestTask"
            };
            task.GetType().GetProperty("Parameters").SetValue(task, parameters);
            var digestionCountDictionary = new Dictionary<(string Accession, string BaseSeqeunce), int>
            {
                { ("Protein1", "SEQUENCE1"), 5 },
                { ("Protein2", "SEQUENCE2"), 10 },
                { ("Protein3", "SEQUENCE3"), 5 }
            };
            task.GetType().GetProperty("DigestionCountDictionary", BindingFlags.NonPublic | BindingFlags.Instance).SetValue(task, digestionCountDictionary);

            // Act
            var method = task.GetType().GetMethod("WriteDigestionCountHistogram", BindingFlags.NonPublic | BindingFlags.Instance);
            method.Invoke(task, null);

            // Assert
            var expectedFilePath = Path.Combine(parameters.OutputFolder, "DigestionCountHistogram.tsv");
            Assert.That(File.Exists(expectedFilePath), Is.True);
            var lines = File.ReadAllLines(expectedFilePath);
            Assert.That(lines.Length, Is.EqualTo(3));
            Assert.That(lines[0], Is.EqualTo("Digestion Products\tCount of Proteins"));
            Assert.That(lines[1], Is.EqualTo("5\t2"));
            Assert.That(lines[2], Is.EqualTo("10\t1"));

            // Cleanup
            Directory.Delete(outputDirectory, true);
        }

        public record DigestionCountTestCase(string DbPath, int MaxIsoforms, bool UseVariableMods, string Name)
        {
            public override string ToString()
            {
                return Name;
            }
        };

        public static IEnumerable<DigestionCountTestCase> GetDigestionCountTestCases()
        {
            // single protein, single peptide
            yield return new DigestionCountTestCase("DatabaseTests//ProteaseModTest.fasta", 1, false, "SingleProteinSinglePeptide_NoMods");
            yield return new DigestionCountTestCase("DatabaseTests//ProteaseModTest.fasta", 1, true, "SingleProteinSinglePeptide_WithMods");
            yield return new DigestionCountTestCase("DatabaseTests//ProteaseModTest.fasta", 128, false, "SingleProteinSinglePeptide_ManyIsoforms_NoMods");
            yield return new DigestionCountTestCase("DatabaseTests//ProteaseModTest.fasta", 128, true, "SingleProteinSinglePeptide_ManyIsoforms_WithMods");

            // single protein, two peptide
            yield return new DigestionCountTestCase("indexEngineTestFasta.fasta", 1, false, "SingleProteinTwoPeptide_NoMods");
            yield return new DigestionCountTestCase("indexEngineTestFasta.fasta", 1, true, "SingleProteinTwoPeptide_WithMods");
            yield return new DigestionCountTestCase("indexEngineTestFasta.fasta", 128, false, "SingleProteinTwoPeptide_ManyIsoforms_NoMods");
            yield return new DigestionCountTestCase("indexEngineTestFasta.fasta", 128, true, "SingleProteinTwoPeptide_ManyIsoforms_WithMods");

            // single protein, many peptides
            yield return new DigestionCountTestCase("DatabaseTests//Q9UHB6.FASTA", 1, false, "SingleProteinManyPeptides_NoMods");
            yield return new DigestionCountTestCase("DatabaseTests//Q9UHB6.FASTA", 1, true, "SingleProteinManyPeptides_WithMods");
            yield return new DigestionCountTestCase("DatabaseTests//Q9UHB6.FASTA", 128, false, "SingleProteinManyPeptides_ManyIsoforms_NoMods");
            yield return new DigestionCountTestCase("DatabaseTests//Q9UHB6.FASTA", 128, true, "SingleProteinManyPeptides_ManyIsoforms_WithMods");

            // many proteins, even more peptides
            yield return new DigestionCountTestCase("TestData//DbForPrunedDb.fasta", 1, false, "ManyProteinsManyPeptides_NoMods");
            yield return new DigestionCountTestCase("TestData//DbForPrunedDb.fasta", 1, true, "ManyProteinsManyPeptides_WithMods");
            yield return new DigestionCountTestCase("TestData//DbForPrunedDb.fasta", 1024, false, "ManyProteinsManyPeptides_ManyIsoforms_NoMods");
            yield return new DigestionCountTestCase("TestData//DbForPrunedDb.fasta", 1024, true, "ManyProteinsManyPeptides_ManyIsoforms_WithMods");
        }

        [Test]
        [TestCaseSource(nameof(GetDigestionCountTestCases))]
        public static void WriteDigestionCountFiles_IsCorrectFromSearchTask(DigestionCountTestCase testCase)
        {
            // Arrange
            string outDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "DigestionCountTest");
            if (Directory.Exists(outDirectory)) 
                Directory.Delete(outDirectory, true);

            var variableMods = testCase.UseVariableMods
                ? new List<(string, string)>
                {
                    ("Common Variable", "Oxidation on M"), ("Common Biological", "Acetylation on A"),
                    ("Common Biological", "Acetylation on G"), ("Common Biological", "Acetylation on K"),
                    ("Common Biological", "Acetylation on M"), ("Common Biological", "Acetylation on P"),
                    ("Common Biological", "Acetylation on S"), ("Common Biological", "Acetylation on T"),
                    ("Common Biological", "Acetylation on X"), ("Common Biological", "Carboxylation on D"),
                    ("Common Biological", "Carboxylation on E"), ("Common Biological", "Carboxylation on K"),
                    ("Common Biological", "Crotonylation on K"), ("Common Biological", "Dimethylation on K"),
                    ("Common Biological", "Dimethylation on R"), ("Common Biological", "Formylation on K"),
                    ("Common Biological", "HexNAc on Nxs"), ("Common Biological", "HexNAc on Nxt"),
                    ("Common Biological", "HexNAc on S"), ("Common Biological", "HexNAc on T"),
                    ("Common Biological", "Hydroxylation on K"), ("Common Biological", "Hydroxylation on N"),
                    ("Common Biological", "Hydroxylation on P"), ("Common Biological", "Methylation on K"),
                    ("Common Biological", "Methylation on Q"), ("Common Biological", "Methylation on R"),
                    ("Common Biological", "Phosphorylation on S"), ("Common Biological", "Phosphorylation on T"),
                    ("Common Biological", "Phosphorylation on Y"), ("Common Biological", "Sulfonation on Y"),
                    ("Common Biological", "Trimethylation on K")
                }
                : [];

            string searchTaskId = "test";
            DigestionParams digestionParams = new DigestionParams(maxModificationIsoforms: testCase.MaxIsoforms, maxMissedCleavages: 0, minPeptideLength: 3);
            var db = new List<DbForTask>() { new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, testCase.DbPath), false) };
            var files = new List<string>() { Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sliced_b6.mzML"), Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "Q9UHB6_Chym_snip.mzML") };
            var tasks = new List<(string, MetaMorpheusTask)>{ (searchTaskId, new SearchTask
            {
                CommonParameters = new CommonParameters(digestionParams: digestionParams, listOfModsVariable: variableMods),
                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    SearchType = SearchType.Classic,
                    SearchTarget = true,
                    DecoyType = DecoyType.None,
                    WriteDigestionProductCountFile = true
                },
            })};

            // convert string modifications to Modification
            object[] parameters = new object[] { "taskId", null, null, null };
            var modConversionMethod = typeof(MetaMorpheusTask).GetMethod("LoadModifications", BindingFlags.NonPublic | BindingFlags.Instance);
            modConversionMethod!.Invoke(tasks.First().Item2, parameters);
            List<Modification> variableModifications = (List<Modification>)parameters[1];

            // Act
            var runner = new EverythingRunnerEngine(tasks, files, db, outDirectory);
            runner.Run();

            // Pull Results from files and calculate from digestion
            var proteins = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, testCase.DbPath), true, DecoyType.None, false, out var errors);
            var digestionResults = proteins.SelectMany(p => p.Digest(digestionParams, [], variableModifications))
                .GroupBy(p => (p.Parent.Accession, p.BaseSequence))
                .ToDictionary(p => p.Key, p => p.ToArray());
            var digestionHistResults = digestionResults.GroupBy(p => p.Value.Length)
                .ToDictionary(p => p.Key, p => p.Count());
            var byProteinLines = File.ReadAllLines(Path.Combine(outDirectory, searchTaskId, "DigestionCountsByProteins.tsv"));
            var histogramLines = File.ReadAllLines(Path.Combine(outDirectory, searchTaskId, "DigestionCountHistogram.tsv"));

            // Assert
            Assert.That(byProteinLines.Length, Is.EqualTo(digestionResults.Count + 1));
            for (int i = 1; i < byProteinLines.Length; i++)
            {
                var split = byProteinLines[i].Split('\t');
                Assert.That(split.Length, Is.EqualTo(3));

                var writtenAccession = split[0];
                var writtenSequence = split[1];
                var writtenCount = int.Parse(split[2]);

                Assert.That(writtenCount, Is.EqualTo(digestionResults[(writtenAccession, writtenSequence)].Length));
            }

            Assert.That(histogramLines.Length, Is.EqualTo(digestionHistResults.Count + 1));
            for (int i = 1; i < histogramLines.Length; i++)
            {
                var split = histogramLines[i].Split('\t');
                Assert.That(split.Length, Is.EqualTo(2));

                var writtenDigestionCount = int.Parse(split[0]);
                var writtenProteinCount = int.Parse(split[1]);

                Assert.That(writtenProteinCount, Is.EqualTo(digestionHistResults[writtenDigestionCount]));
            }

            // Cleanup
            Directory.Delete(outDirectory, true);
        }
        [Test]
        public static void WriteDigestionCountFiles_DoesNotIncludeDecoys_WhenNotIntended()
        {
            // Arrange
            var task = new PostSearchAnalysisTask();
            var outputDirectory = Path.Combine(TestContext.CurrentContext.WorkDirectory, "DigestionCountTest");
            if (Directory.Exists(outputDirectory))
                Directory.Delete(outputDirectory, true);
            Directory.CreateDirectory(outputDirectory);
            var parameters = new PostSearchAnalysisParameters
            {
                OutputFolder = outputDirectory,
                SearchTaskId = "TestTask",
                SearchParameters = new SearchParameters
                {
                    WriteDecoys = false
                }
            };

            task.GetType().GetProperty("Parameters").SetValue(task, parameters);
            var digestionCountDictionary = new Dictionary<(string Accession, string BaseSeqeunce), int>
            {
                { ("DECOY_Protein1", "SEQUENCE1"), 5 },
                { ("Protein2", "SEQUENCE2"), 10 }
            };
            task.GetType().GetProperty("DigestionCountDictionary", BindingFlags.NonPublic | BindingFlags.Instance).SetValue(task, digestionCountDictionary);

            // Act
            var method = task.GetType().GetMethod("WriteDigestionCountByProtein", BindingFlags.NonPublic | BindingFlags.Instance);
            method!.Invoke(task, null);

            // Assert
            var expectedFilePath = Path.Combine(parameters.OutputFolder, "DigestionCountsByProteins.tsv");
            Assert.That(File.Exists(expectedFilePath), Is.True);

            var lines = File.ReadAllLines(expectedFilePath);
            Assert.That(lines.Length, Is.EqualTo(2));
            Assert.That(lines[0], Is.EqualTo("Protein Accession\tPrimary Sequence\tDigestion Products"));
            Assert.That(lines[1], Is.EqualTo("Protein2\tSEQUENCE2\t10"));

            // Cleanup
            Directory.Delete(parameters.OutputFolder, true);
        }

        [Test]
        public static void WriteDigestionCountFiles_IncludesDecoys_WhenIntended()
        {
            // Arrange
            var task = new PostSearchAnalysisTask();
            var outputDirectory = Path.Combine(TestContext.CurrentContext.WorkDirectory, "DigestionCountTest");
            if (Directory.Exists(outputDirectory))
                Directory.Delete(outputDirectory, true);
            Directory.CreateDirectory(outputDirectory);
            var parameters = new PostSearchAnalysisParameters
            {
                OutputFolder = outputDirectory,
                SearchTaskId = "TestTask",
                SearchParameters = new SearchParameters
                {
                    WriteDecoys = true
                }
            };

            task.GetType().GetProperty("Parameters").SetValue(task, parameters);
            var digestionCountDictionary = new Dictionary<(string Accession, string BaseSeqeunce), int>
            {
                { ("DECOY_Protein1", "SEQUENCE1"), 5 },
                { ("Protein2", "SEQUENCE2"), 10 }
            };
            task.GetType().GetProperty("DigestionCountDictionary", BindingFlags.NonPublic | BindingFlags.Instance).SetValue(task, digestionCountDictionary);

            // Act
            var method = task.GetType().GetMethod("WriteDigestionCountByProtein", BindingFlags.NonPublic | BindingFlags.Instance);
            method!.Invoke(task, null);

            // Assert
            var expectedFilePath = Path.Combine(parameters.OutputFolder, "DigestionCountsByProteins.tsv");
            Assert.That(File.Exists(expectedFilePath), Is.True);

            var lines = File.ReadAllLines(expectedFilePath);
            Assert.That(lines.Length, Is.EqualTo(3));
            Assert.That(lines[0], Is.EqualTo("Protein Accession\tPrimary Sequence\tDigestion Products"));
            Assert.That(lines[1], Is.EqualTo("DECOY_Protein1\tSEQUENCE1\t5"));
            Assert.That(lines[2], Is.EqualTo("Protein2\tSEQUENCE2\t10"));

            // Cleanup
            Directory.Delete(parameters.OutputFolder, true);
        }

        /// <summary>
        /// A resolved sequence does not imply a resolved mass. Two best matches whose modifications
        /// share an IdWithMotif but not a mass stringify to one FullSequence, so the PSM looks
        /// unambiguous to every check FilteredPsms.Filter makes, while its mass stays null.
        /// </summary>
        [Test]
        public static void ModMassCollisionLeavesMassUnresolvedButPassesTheQuantFilter()
        {
            CommonParameters commonParameters = new();
            var (light, heavy) = TwoPeptidesSharingAFullSequence();

            Assert.That(heavy.FullSequence, Is.EqualTo(light.FullSequence));
            Assert.That(heavy.MonoisotopicMass, Is.Not.EqualTo(light.MonoisotopicMass).Within(1e-9));

            SpectralMatch psm = new PeptideSpectralMatch(light, 0, 10, 0,
                NullMassTestScan("fake.mzML", commonParameters, 1), commonParameters, new List<MatchedFragmentIon>());
            psm.AddOrReplace(heavy, 10, 0, reportAllAmbiguity: true, new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            psm.SetFdrValues(1, 0, 0.0, 1, 0, 0.0, 0, 0.0);

            Assert.That(psm.BaseSequence, Is.Not.Null.And.Not.Empty);
            Assert.That(psm.FullSequence, Is.Not.Null.And.Not.Empty);
            Assert.That(psm.BioPolymerWithSetModsMonoisotopicMass, Is.Null);

            var survivors = FilteredPsms.Filter(new List<SpectralMatch> { psm }, commonParameters,
                includeDecoys: false, includeContaminants: true, includeAmbiguous: false,
                includeAmbiguousMods: false, includeHighQValuePsms: false);

            Assert.That(survivors.Count(), Is.EqualTo(1));
        }

        /// <summary>
        /// Residue.ResidueMonoisotopicMass is double.NaN for B, X and Z. NaN fails the tolerance
        /// comparison in PsmTsvWriter.Resolve, so even a single unambiguous match resolves to a null mass.
        /// </summary>
        [Test]
        public static void AmbiguousResidueLeavesMassUnresolved()
        {
            CommonParameters commonParameters = new();
            var protein = new Protein("PEPTXDEK", "ACC_X");
            var digestionParams = new DigestionParams(protease: "trypsin", minPeptideLength: 1);
            var peptide = (PeptideWithSetModifications)protein
                .Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();

            Assert.That(double.IsNaN(peptide.MonoisotopicMass), Is.True);

            SpectralMatch psm = new PeptideSpectralMatch(peptide, 0, 10, 0,
                NullMassTestScan("fake.mzML", commonParameters, 1), commonParameters, new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();

            Assert.That(psm.BestMatchingBioPolymersWithSetMods.Count(), Is.EqualTo(1));
            Assert.That(psm.FullSequence, Is.Not.Null.And.Not.Empty);
            Assert.That(psm.BioPolymerWithSetModsMonoisotopicMass, Is.Null);
        }

        /// <summary>
        /// RemovePsmsWithoutResolvedMass is public API and was previously reachable only through
        /// QuantificationAnalysis by reflection. Both halves of its contract are asserted here: the
        /// null-mass PSM goes, the resolved one stays, and the return value is the number dropped.
        /// </summary>
        [Test]
        public static void RemovePsmsWithoutResolvedMass_DropsOnlyTheMasslessOnes()
        {
            CommonParameters commonParameters = new();

            var protein = new Protein("PEPTIDEKPEPTIDER", "ACC_DIRECT");
            var digestionParams = new DigestionParams(protease: "trypsin", minPeptideLength: 1);
            var resolved = (PeptideWithSetModifications)protein
                .Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .First(p => p.BaseSequence == "PEPTIDER");

            SpectralMatch withMass = new PeptideSpectralMatch(resolved, 0, 10, 0,
                NullMassTestScan("fake.mzML", commonParameters, 1), commonParameters, new List<MatchedFragmentIon>());
            withMass.ResolveAllAmbiguities();
            withMass.SetFdrValues(1, 0, 0.0, 1, 0, 0.0, 0, 0.0);

            var (light, heavy) = TwoPeptidesSharingAFullSequence();
            SpectralMatch withoutMass = new PeptideSpectralMatch(light, 0, 10, 1,
                NullMassTestScan("fake.mzML", commonParameters, 2), commonParameters, new List<MatchedFragmentIon>());
            withoutMass.AddOrReplace(heavy, 10, 0, reportAllAmbiguity: true, new List<MatchedFragmentIon>());
            withoutMass.ResolveAllAmbiguities();
            withoutMass.SetFdrValues(1, 0, 0.0, 1, 0, 0.0, 0, 0.0);

            Assert.That(withMass.BioPolymerWithSetModsMonoisotopicMass, Is.Not.Null);
            Assert.That(withoutMass.BioPolymerWithSetModsMonoisotopicMass, Is.Null);

            var filtered = FilteredPsms.Filter(new List<SpectralMatch> { withMass, withoutMass }, commonParameters,
                includeDecoys: false, includeContaminants: true, includeAmbiguous: false,
                includeAmbiguousMods: false, includeHighQValuePsms: false);
            Assert.That(filtered.Count(), Is.EqualTo(2), "both PSMs must survive sequence filtering, or this proves nothing");

            int dropped = filtered.RemovePsmsWithoutResolvedMass();

            Assert.Multiple(() =>
            {
                Assert.That(dropped, Is.EqualTo(1));
                Assert.That(filtered.Count(), Is.EqualTo(1));
                Assert.That(filtered.Single(), Is.SameAs(withMass), "the PSM that has a mass must be the one kept");
            });
        }

        /// <summary>
        /// The other half of the guard's stated behaviour: a run where every mass resolves drops
        /// nothing and emits no warning. Every other new test feeds at least one null-mass PSM, so an
        /// inverted guard, or one returning a nonzero count when it removed nothing, would pass them all.
        /// </summary>
        [Test]
        public static void RemovePsmsWithoutResolvedMass_KeepsEverythingWhenAllMassesResolve()
        {
            CommonParameters commonParameters = new();

            var protein = new Protein("PEPTIDEKPEPTIDER", "ACC_ALLGOOD");
            var digestionParams = new DigestionParams(protease: "trypsin", minPeptideLength: 1);
            var peptides = protein.Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .Cast<PeptideWithSetModifications>().ToList();

            var psms = new List<SpectralMatch>();
            int scan = 1;
            foreach (var peptide in peptides)
            {
                SpectralMatch psm = new PeptideSpectralMatch(peptide, 0, 10, scan - 1,
                    NullMassTestScan("fake.mzML", commonParameters, scan++), commonParameters, new List<MatchedFragmentIon>());
                psm.ResolveAllAmbiguities();
                psm.SetFdrValues(1, 0, 0.0, 1, 0, 0.0, 0, 0.0);
                Assert.That(psm.BioPolymerWithSetModsMonoisotopicMass, Is.Not.Null);
                psms.Add(psm);
            }

            var filtered = FilteredPsms.Filter(psms, commonParameters,
                includeDecoys: false, includeContaminants: true, includeAmbiguous: false,
                includeAmbiguousMods: false, includeHighQValuePsms: false);
            int before = filtered.Count();

            int dropped = filtered.RemovePsmsWithoutResolvedMass();

            Assert.Multiple(() =>
            {
                Assert.That(dropped, Is.Zero, "nothing should be dropped when every mass resolved");
                Assert.That(filtered.Count(), Is.EqualTo(before));
            });
        }

        /// <summary>
        /// One PSM without a resolvable mass used to throw out of QuantificationAnalysis, which
        /// EngineCrashed turns into a Quantification_crash.txt and a null FlashLfqResults -- losing
        /// quantification for the whole run. It should cost only the offending PSM, plus a warning.
        /// </summary>
        [Test]
        public static void QuantificationSkipsPsmsWithUnresolvedMass()
        {
            CommonParameters commonParameters = new();
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "NullMassQuant");
            if (Directory.Exists(outputFolder))
            {
                Directory.Delete(outputFolder, true);
            }
            Directory.CreateDirectory(outputFolder);

            string mzmlPath = Path.Combine(outputFolder, "fake.mzML");
            var protein = new Protein("PEPTIDEKPEPTIDERPEPTIDEK", "ACC_QUANT");
            var digestionParams = new DigestionParams(protease: "trypsin", minPeptideLength: 1);
            var digestionProducts = protein.Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .Cast<PeptideWithSetModifications>().ToList();

            SpectralMatch quantifiablePsm = new PeptideSpectralMatch(
                digestionProducts.First(p => p.BaseSequence == "PEPTIDER"), 0, 10, 0,
                NullMassTestScan(mzmlPath, commonParameters, 1), commonParameters, new List<MatchedFragmentIon>());
            quantifiablePsm.ResolveAllAmbiguities();
            quantifiablePsm.SetFdrValues(1, 0, 0.0, 1, 0, 0.0, 0, 0.0);

            var (light, heavy) = TwoPeptidesSharingAFullSequence();
            SpectralMatch masslessPsm = new PeptideSpectralMatch(light, 0, 10, 1,
                NullMassTestScan(mzmlPath, commonParameters, 2), commonParameters, new List<MatchedFragmentIon>());
            masslessPsm.AddOrReplace(heavy, 10, 0, reportAllAmbiguity: true, new List<MatchedFragmentIon>());
            masslessPsm.ResolveAllAmbiguities();
            masslessPsm.SetFdrValues(1, 0, 0.0, 1, 0, 0.0, 0, 0.0);

            Assert.That(quantifiablePsm.BioPolymerWithSetModsMonoisotopicMass, Is.Not.Null);
            Assert.That(masslessPsm.BioPolymerWithSetModsMonoisotopicMass, Is.Null);

            // The MS1 signal has to correspond to the surviving identification, or FlashLFQ finds
            // nothing and every intensity comes back zero -- which would let this test pass while
            // proving only that an Identification was constructed.
            WriteMs1FixtureFor(mzmlPath, "PEPTIDER");

            List<string> warnings = new();
            void Handler(object sender, StringEventArgs e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += Handler;

            try
            {
                PostSearchAnalysisParameters parameters = new()
                {
                    SearchParameters = new SearchParameters
                    {
                        DoLabelFreeQuantification = true,
                        DoMultiplexQuantification = false,
                        MatchBetweenRuns = false,
                        Normalize = false,
                    },
                    OutputFolder = outputFolder,
                    IndividualResultsOutputFolder = outputFolder,
                    SearchTaskId = "TestTask",
                    AllSpectralMatches = new List<SpectralMatch> { quantifiablePsm, masslessPsm },
                    CurrentRawFileList = new List<string> { mzmlPath },
                    MyFileManager = new MyFileManager(true),
                    FixedModifications = new List<Modification>(),
                    VariableModifications = new List<Modification>(),
                    ListOfDigestionParams = new HashSet<IDigestionParams> { digestionParams },
                    DatabaseFilenameList = new List<DbForTask>(),
                    FileSettingsList = new FileSpecificParameters[] { null },
                };

                PostSearchAnalysisTask task = new()
                {
                    Parameters = parameters,
                    CommonParameters = commonParameters,
                    FileSpecificParameters = new List<(string, CommonParameters)> { (mzmlPath, commonParameters) },
                };

                typeof(PostSearchAnalysisTask)
                    .GetMethod("QuantificationAnalysis", BindingFlags.NonPublic | BindingFlags.Instance)!
                    .Invoke(task, null);

                Assert.That(File.Exists(Path.Combine(outputFolder, "Quantification_crash.txt")), Is.False);
                // The count, not just the tail of the message: exactly one PSM lacked a mass, and a
                // method that dropped the wrong number would still match on the sentence alone.
                Assert.That(warnings.Any(w => w.Contains("1 PSM(s) were excluded from quantification")), Is.True);
                Assert.That(parameters.FlashLfqResults, Is.Not.Null);
                Assert.That(parameters.FlashLfqResults.PeptideModifiedSequences.Count, Is.EqualTo(1));

                // ...and that the survivor was actually quantified, rather than merely counted. Without
                // a matching MS1 peak this is zero while every other assertion above still passes.
                var survivor = parameters.FlashLfqResults.PeptideModifiedSequences.Single().Value;
                var quantFile = parameters.FlashLfqResults.SpectraFiles.Single();
                Assert.That(survivor.GetIntensity(quantFile), Is.GreaterThan(0),
                    "the surviving PSM must be quantified, not just present in the results");
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= Handler;
                Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// The SILAC path rebuilds the quantification list after the first guard has already run:
        /// SetSilacFilteredPsms replaces it wholesale with clones synthesised from labeled base
        /// sequences, so a PSM whose mass resolved before the block can be replaced by one whose mass
        /// does not. The label residue here has no defined mass -- the same NaN route as B, X and Z --
        /// which is the state PostSearchAnalysisTask sees whenever the label residues were never added
        /// to Residue's process-wide table, since it is SearchTask, not this task, that adds them.
        /// Without the second guard that clone reaches BioPolymerWithSetModsMonoisotopicMass.Value and
        /// takes quantification for the whole run down with it.
        /// </summary>
        [Test]
        public static void SilacQuantificationSkipsLabeledPsmsWithUnresolvedMass()
        {
            // 'X' is one of GlobalVariables.InvalidAminoAcids, so SilacConversions.UpdateAminoAcidLabel
            // never assigns it and nothing in the suite registers a mass for it. The label residue is
            // therefore massless no matter what other tests have added to Residue's static table, which
            // keeps this test independent of run order.
            const char masslessLabelResidue = 'X';
            Assert.That(double.IsNaN(Proteomics.AminoAcidPolymer.Residue.ResidueMonoisotopicMass[masslessLabelResidue]), Is.True,
                "the label residue must have no defined mass, or the labeled clone would resolve and this test would prove nothing");

            CommonParameters commonParameters = new();
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "SilacNullMassQuant");
            if (Directory.Exists(outputFolder))
            {
                Directory.Delete(outputFolder, true);
            }
            Directory.CreateDirectory(outputFolder);

            string mzmlPath = Path.Combine(outputFolder, "fake.mzML");
            var protein = new Protein("PEPTIDEKPEPTIDER", "ACC_SILAC");

            // GeneratehUnlabeledProteinsForSilac is what puts the light form into the SILAC list
            // alongside the labeled clone, so the run has something left to quantify once the clone is
            // dropped -- which is the behaviour under test, not merely that nothing threw.
            var digestionParams = new DigestionParams(protease: "trypsin", minPeptideLength: 1,
                generateUnlabeledProteinsForSilac: true);
            var digestionProducts = protein.Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .Cast<PeptideWithSetModifications>().ToList();

            // PEPTIDEK carries the labeled residue, so it gets a labeled clone. PEPTIDER does not, so
            // SILAC leaves it as the light form only: one PSM must be dropped, not both.
            SpectralMatch lysinePsm = ResolvedPsm(digestionProducts.First(p => p.BaseSequence == "PEPTIDEK"),
                mzmlPath, commonParameters, 1);
            SpectralMatch argininePsm = ResolvedPsm(digestionProducts.First(p => p.BaseSequence == "PEPTIDER"),
                mzmlPath, commonParameters, 2);

            // Both masses resolve, so the guard ahead of the SILAC block has nothing to do and only the
            // second guard can account for what this test observes.
            Assert.That(lysinePsm.BioPolymerWithSetModsMonoisotopicMass, Is.Not.Null);
            Assert.That(argininePsm.BioPolymerWithSetModsMonoisotopicMass, Is.Not.Null);

            WriteMs1FixtureFor(mzmlPath, "PEPTIDEK", "PEPTIDER");

            List<string> warnings = new();
            void Handler(object sender, StringEventArgs e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += Handler;

            try
            {
                PostSearchAnalysisParameters parameters = new()
                {
                    SearchParameters = new SearchParameters
                    {
                        DoLabelFreeQuantification = true,
                        DoMultiplexQuantification = false,
                        MatchBetweenRuns = false,
                        Normalize = false,
                        SilacLabels = new List<SilacLabel>
                        {
                            new SilacLabel('K', masslessLabelResidue, "C{13}6", 6.020129)
                        },
                    },
                    OutputFolder = outputFolder,
                    IndividualResultsOutputFolder = outputFolder,
                    SearchTaskId = "TestTask",
                    AllSpectralMatches = new List<SpectralMatch> { lysinePsm, argininePsm },
                    CurrentRawFileList = new List<string> { mzmlPath },
                    MyFileManager = new MyFileManager(true),
                    FixedModifications = new List<Modification>(),
                    VariableModifications = new List<Modification>(),
                    ListOfDigestionParams = new HashSet<IDigestionParams> { digestionParams },
                    DatabaseFilenameList = new List<DbForTask>(),
                    FileSettingsList = new FileSpecificParameters[] { null },
                };

                PostSearchAnalysisTask task = new()
                {
                    Parameters = parameters,
                    CommonParameters = commonParameters,
                    FileSpecificParameters = new List<(string, CommonParameters)> { (mzmlPath, commonParameters) },
                };

                typeof(PostSearchAnalysisTask)
                    .GetMethod("QuantificationAnalysis", BindingFlags.NonPublic | BindingFlags.Instance)!
                    .Invoke(task, null);

                Assert.Multiple(() =>
                {
                    Assert.That(File.Exists(Path.Combine(outputFolder, "Quantification_crash.txt")), Is.False,
                        "the labeled clone must not take the whole quantification down with it");

                    // The SILAC message, and the count in it: exactly one clone lacked a mass.
                    Assert.That(warnings.Any(w => w.Contains("1 SILAC PSM(s) were excluded from quantification")), Is.True);

                    // ...and not the message from the guard ahead of the SILAC block, which had nothing
                    // to drop. Asserting on the shared tail alone would let either guard satisfy this.
                    Assert.That(warnings.Any(w => w.Contains("PSM(s) were excluded") && !w.Contains("SILAC")), Is.False,
                        "the first guard saw two resolved masses and must not have dropped anything");

                    Assert.That(parameters.FlashLfqResults, Is.Not.Null);
                });

                // The labeled sequence is gone and both light forms are still quantified: dropping one
                // clone cost that clone, not the run.
                var quantifiedSequences = parameters.FlashLfqResults.PeptideModifiedSequences.Keys.ToList();
                Assert.Multiple(() =>
                {
                    Assert.That(quantifiedSequences, Does.Not.Contain("PEPTIDE" + masslessLabelResidue));
                    Assert.That(quantifiedSequences, Does.Contain("PEPTIDEK"));
                    Assert.That(quantifiedSequences, Does.Contain("PEPTIDER"));
                });

                // Present in the results is not the same as measured. Without a matching MS1 envelope
                // every intensity is zero while every assertion above still passes.
                foreach (string sequence in new[] { "PEPTIDEK", "PEPTIDER" })
                {
                    var peptide = parameters.FlashLfqResults.PeptideModifiedSequences[sequence];
                    Assert.That(parameters.FlashLfqResults.SpectraFiles.Any(f => peptide.GetIntensity(f) > 0), Is.True,
                        sequence + " survived the guard but was never quantified");
                }
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= Handler;
                Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// A PSM on a single unambiguous peptide, resolved and given passing FDR values, so it reaches
        /// quantification with a mass.
        /// </summary>
        private static SpectralMatch ResolvedPsm(PeptideWithSetModifications peptide, string mzmlPath,
            CommonParameters commonParameters, int scanNumber)
        {
            SpectralMatch psm = new PeptideSpectralMatch(peptide, 0, 10, scanNumber - 1,
                NullMassTestScan(mzmlPath, commonParameters, scanNumber), commonParameters, new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            psm.SetFdrValues(1, 0, 0.0, 1, 0, 0.0, 0, 0.0);
            return psm;
        }

        /// <summary>
        /// Writes an mzML whose MS1 scans carry the real isotopic envelope of each named peptide at
        /// charge 2. A single monoisotopic peak is not enough -- FlashLFQ accepts a feature only when it
        /// can build an envelope -- and the envelopes are derived from the sequences rather than
        /// hardcoded, so the fixture cannot drift away from the identifications the tests build. Three
        /// scans span the identifications' retention time (10.0) so the peak has a shape to integrate
        /// across rather than a single point.
        /// </summary>
        private static void WriteMs1FixtureFor(string mzmlPath, params string[] baseSequences)
        {
            var envelopes = baseSequences
                .Select(sequence =>
                {
                    ChemicalFormula formula = new Proteomics.AminoAcidPolymer.Peptide(sequence).GetChemicalFormula();
                    return IsotopicDistribution.GetDistribution(formula, 0.125, 1e-8);
                })
                .ToList();

            var scans = new[] { 9.8, 10.0, 10.2 }
                .Select((rt, i) =>
                {
                    double scale = i == 1 ? 1e6 : 5e5;
                    var peaks = envelopes
                        .SelectMany(e => e.Masses.Select(m => m.ToMz(2)).Zip(e.Intensities, (mz, intensity) => (mz, intensity: intensity * scale)))
                        .OrderBy(peak => peak.mz)
                        .ToList();
                    double[] mzs = peaks.Select(peak => peak.mz).ToArray();
                    double[] intensities = peaks.Select(peak => peak.intensity).ToArray();
                    return new MsDataScan(
                        new MzSpectrum(mzs, intensities, false),
                        i + 1, 1, true, Polarity.Positive, rt, new MzRange(0, 2000), "f",
                        MZAnalyzerType.Orbitrap, intensities.Sum(), 1.0, null, "scan=" + (i + 1));
                })
                .ToArray();

            var dataFile = new GenericMsDataFile(scans, new SourceFile(null, null, null, null, null));
            Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(dataFile, mzmlPath, false);
        }

        private static Ms2ScanWithSpecificMass NullMassTestScan(string filePath, CommonParameters commonParameters, int scanNumber)
        {
            MsDataScan scan = new(
                new MzSpectrum(new[] { 100.0, 200.0 }, new[] { 1.0, 1.0 }, false),
                oneBasedScanNumber: scanNumber, msnOrder: 2, isCentroid: true,
                polarity: Polarity.Positive, retentionTime: 10.0,
                scanWindowRange: new MzRange(0, 2000), scanFilter: "f",
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 2,
                injectionTime: 1.0, noiseData: null, nativeId: $"scan={scanNumber}",
                selectedIonMz: 500.0, selectedIonChargeStateGuess: 2,
                selectedIonIntensity: 1, isolationMZ: 500.0, isolationWidth: 2,
                dissociationType: DissociationType.HCD, oneBasedPrecursorScanNumber: null,
                selectedIonMonoisotopicGuessMz: 500.0);

            return new Ms2ScanWithSpecificMass(scan, 500.0, 2, filePath, commonParameters);
        }

        /// <summary>
        /// Two peptides differing only in the mass of a modification whose ModificationType and
        /// IdWithMotif are identical. FullSequence is built from those two strings alone.
        /// </summary>
        private static (PeptideWithSetModifications Light, PeptideWithSetModifications Heavy) TwoPeptidesSharingAFullSequence()
        {
            ModificationMotif.TryGetMotif("P", out var motif);
            Modification lightMod = new(_originalId: "Collide", _modificationType: "TestType",
                _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 10.0);
            Modification heavyMod = new(_originalId: "Collide", _modificationType: "TestType",
                _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 20.0);

            var protein = new Protein("PEPTIDEK", "ACC_COLLIDE");
            var digestionParams = new DigestionParams(protease: "trypsin", minPeptideLength: 1);
            var unmodified = (PeptideWithSetModifications)protein
                .Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();

            PeptideWithSetModifications WithMod(Modification mod) => new(protein, digestionParams,
                unmodified.OneBasedStartResidue, unmodified.OneBasedEndResidue,
                unmodified.CleavageSpecificityForFdrCategory, unmodified.PeptideDescription,
                unmodified.MissedCleavages, new Dictionary<int, Modification> { { 2, mod } },
                unmodified.NumFixedMods);

            return (WithMod(lightMod), WithMod(heavyMod));
        }

    }
}