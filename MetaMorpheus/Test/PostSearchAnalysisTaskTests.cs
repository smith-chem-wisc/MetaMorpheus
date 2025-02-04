using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Reflection;
using EngineLayer;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
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
    }
}