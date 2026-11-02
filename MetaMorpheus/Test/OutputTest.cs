using Nett;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using EngineLayer.DatabaseLoading;
using MetaMorpheusCommandLine;
using NUnit.Framework.Legacy;
using Readers;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class OutputTest
    {
        [Test]
        public static void TestIndividualFileOutput()
        {
            string subFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"IndividualOutputTest");
            if(Directory.Exists(subFolder)) 
                Directory.Delete(subFolder, true);
            Directory.CreateDirectory(subFolder);
            string outputFolder = Path.Combine(subFolder, "Results");
            SearchTask allowFilesTask = new SearchTask();
            allowFilesTask.SearchParameters.WriteIndividualFiles = true;
            allowFilesTask.SearchParameters.CompressIndividualFiles = false;
            allowFilesTask.SearchParameters.WriteMzId = true;

            SearchTask compressFilesTask = new SearchTask();
            compressFilesTask.SearchParameters.WriteIndividualFiles = true;
            compressFilesTask.SearchParameters.CompressIndividualFiles = true;
            compressFilesTask.SearchParameters.WriteMzId = true;

            SearchTask noFilesTask = new SearchTask();
            noFilesTask.SearchParameters.WriteIndividualFiles = false;
            noFilesTask.SearchParameters.WriteMzId = true;

            PeptideWithSetModifications pwsm = new PeptideWithSetModifications("AAFNSGK", null);
            List<(string, MetaMorpheusTask)> tasks = new List<(string, MetaMorpheusTask)> { ("allowFiles", allowFilesTask), ("compressFiles", compressFilesTask), ("noFiles", noFilesTask) };
            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "gapdh.fasta"), false);

            TestDataFile datafile = new TestDataFile(pwsm);
            string pathOne = Path.Combine(subFolder, "fileOne.mzml");
            string pathTwo = Path.Combine(subFolder, "fileTwo.mzml");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(datafile, pathOne, false);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(datafile, pathTwo, false);

            new EverythingRunnerEngine(tasks, new List<string> { pathOne, pathTwo }, new List<DbForTask> { db }, outputFolder).Run();

            //check that the first task wrote everything fine
            HashSet<string> expectedFiles = new HashSet<string>
            {
                ".mzID",
                "_Peptides.psmtsv",
                "_ProteinGroups.tsv",
                "_PSMs.psmtsv",
                "_PSMsFormattedForPercolator.tab",
                "_QuantifiedPeaks.tsv"
            };
            HashSet<string> writtenFiles = new HashSet<string>(Directory.GetFiles(Path.Combine(outputFolder, "allowFiles", "Individual File Results")).Select(v => Path.GetFileName(v).Substring(7)));
            //check they're the same
            Assert.That(expectedFiles.Except(writtenFiles).Count() == 0);

            //check the second one is compressed and contains all the information
            writtenFiles = new HashSet<string>(Directory.GetFiles(Path.Combine(outputFolder, "compressFiles")).Select(v => Path.GetFileName(v)));
            //check the zip exists
            Assert.That(writtenFiles.Contains("Individual File Results.zip"));
            //check the original folder does not exist
            string[] subfolders = Directory.GetDirectories(Path.Combine(outputFolder, "compressFiles"));
            Assert.That(subfolders.Length == 0);
            ZipFile.ExtractToDirectory(Path.Combine(outputFolder, "compressFiles", "Individual File Results.zip"), Path.Combine(outputFolder, "compressFiles", "Individual File Results"));
            //read the extracted files
            writtenFiles = new HashSet<string>(Directory.GetFiles(Path.Combine(outputFolder, "compressFiles", "Individual File Results")).Select(v => Path.GetFileName(v).Substring(7)));
            //check they're the same
            Assert.That(expectedFiles.Except(writtenFiles).Count() == 0);

            //check the last one to make sure nothing was written except for the mzID files
            writtenFiles = new HashSet<string>(Directory.GetFiles(Path.Combine(outputFolder, "noFiles", "Individual File Results")).Select(v => Path.GetFileName(v).Substring(7)));
            Assert.That(writtenFiles.Count == 1);
            Assert.That(writtenFiles.Contains(".mzID"));

            Directory.Delete(outputFolder, true);

            //Do a check that we don't crash if there's only one file but somebody tries to zip the individual file results
            SearchTask weirdTask = new SearchTask();
            weirdTask.SearchParameters.CompressIndividualFiles = true;
            weirdTask.SearchParameters.WriteMzId = false;
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("weird", weirdTask) }, new List<string> { pathOne }, new List<DbForTask> { db }, outputFolder).Run();
            //check that a zip was not created
            writtenFiles = new HashSet<string>(Directory.GetFiles(Path.Combine(outputFolder, "weird")));
            Assert.That(!writtenFiles.Contains("Individual File Results.zip"));
            Directory.Delete(subFolder, true);
        }

        /// <summary>
        /// WriteIndividualFiles is off by default. A default search of two files writes no per-file PSM, peptide,
        /// protein group or quantified peak table and no per-file counts in results.txt, but the per-file .mzID
        /// still lands in "Individual File Results", because WriteMzId governs it. The default task TOML that
        /// CMD -g writes says false. A task TOML without the key reads as off; one that sets it to true still
        /// writes the per-file tables.
        /// </summary>
        [Test]
        public static void IndividualFileTablesAreOffByDefault()
        {
            Assert.That(new SearchParameters().WriteIndividualFiles, Is.False);

            string subFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"IndividualFilesOffByDefaultTest");
            if (Directory.Exists(subFolder))
                Directory.Delete(subFolder, true);
            Directory.CreateDirectory(subFolder);
            try
            {
                string defaultTomls = Path.Combine(subFolder, "defaultTomls");
                CommandLineSettings.GenerateDefaultTaskTomls(defaultTomls);
                string[] savedLines = File.ReadAllLines(Path.Combine(defaultTomls, "SearchTask.toml"));
                Func<string, bool> isTheKey = line => line.Trim().StartsWith("WriteIndividualFiles ");
                Assert.That(savedLines.Where(isTheKey).Select(line => line.Trim()), Is.EqualTo(new[] { "WriteIndividualFiles = false" }),
                    "the default task TOML carries the key, set to false");

                string keyAbsentToml = Path.Combine(subFolder, "keyAbsent.toml");
                File.WriteAllLines(keyAbsentToml, savedLines.Where(line => !isTheKey(line)));
                string keyTrueToml = Path.Combine(subFolder, "keyTrue.toml");
                File.WriteAllLines(keyTrueToml, savedLines.Select(line => isTheKey(line) ? "WriteIndividualFiles = true" : line));

                SearchTask defaultTask = Toml.ReadFile<SearchTask>(keyAbsentToml, MetaMorpheusTask.tomlConfig);
                SearchTask optedInTask = Toml.ReadFile<SearchTask>(keyTrueToml, MetaMorpheusTask.tomlConfig);
                Assert.That(defaultTask.SearchParameters.WriteIndividualFiles, Is.False);
                Assert.That(optedInTask.SearchParameters.WriteIndividualFiles, Is.True);
                Assert.That(defaultTask.SearchParameters.WriteMzId, Is.True, "premise: .mzID is written by default");

                DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "gapdh.fasta"), false);
                TestDataFile datafile = new TestDataFile(new PeptideWithSetModifications("AAFNSGK", null));
                string pathOne = Path.Combine(subFolder, "fileOne.mzml");
                string pathTwo = Path.Combine(subFolder, "fileTwo.mzml");
                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(datafile, pathOne, false);
                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(datafile, pathTwo, false);

                string outputFolder = Path.Combine(subFolder, "Results");
                var tasks = new List<(string, MetaMorpheusTask)> { ("default", defaultTask), ("optedIn", optedInTask) };
                new EverythingRunnerEngine(tasks, new List<string> { pathOne, pathTwo }, new List<DbForTask> { db }, outputFolder).Run();

                Func<string, string[]> individualFiles = task =>
                    Directory.GetFiles(Path.Combine(outputFolder, task, "Individual File Results")).Select(Path.GetFileName).ToArray();
                Func<string, string[]> perFileCounts = task =>
                    File.ReadAllLines(Path.Combine(outputFolder, task, "results.txt"))
                        .Where(line => line.StartsWith("fileOne - ") || line.StartsWith("fileTwo - ")).ToArray();

                Assert.That(individualFiles("default"), Is.EquivalentTo(new[] { "fileOne.mzID", "fileTwo.mzID" }));
                Assert.That(perFileCounts("default"), Is.Empty);

                Assert.That(individualFiles("optedIn"), Is.SupersetOf(new[]
                {
                    "fileOne.mzID", "fileOne_PSMs.psmtsv", "fileOne_Peptides.psmtsv", "fileOne_ProteinGroups.tsv", "fileOne_QuantifiedPeaks.tsv",
                    "fileTwo.mzID", "fileTwo_PSMs.psmtsv", "fileTwo_Peptides.psmtsv", "fileTwo_ProteinGroups.tsv", "fileTwo_QuantifiedPeaks.tsv"
                }));
                Assert.That(perFileCounts("optedIn"), Has.Some.StartsWith("fileOne - Target PSMs"));
            }
            finally
            {
                Directory.Delete(subFolder, true);
            }
        }

        [Test]
        public static void TestOpairFileOutput()
        {
            var task = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);

            string currentDirectory = Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData");

            Directory.CreateDirectory(currentDirectory);
            DbForTask db = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P16150.fasta"), false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) }, new List<string> { spectraFile }, new List<DbForTask> { db }, currentDirectory).Run();

            string outputDirectory = Path.Combine(currentDirectory, @"Task");

            HashSet<string> writtenFiles = new(Directory.GetFiles(outputDirectory).Select(v => Path.GetFileName(v)));

            //check that the first task wrote everything fine
            HashSet<string> expectedFiles = new()
            {
                "AllPSMs.psmtsv",
                "_AllProteinGroups.tsv",
                "oglyco.psmtsv",
                "protein_oglyco_localization.tsv",
                "seen_oglyco_localization.tsv",
                "results.txt",
                "AutoGeneratedManuscriptProse.txt"
            };
            
            //check they're the same
            CollectionAssert.AreEquivalent(expectedFiles, writtenFiles);    

            Directory.Delete(currentDirectory, true);
        }
    }
}