using EngineLayer;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class OutputTest
    {
        [Test]
        public static void TestIndividualFileOutput()
        {
            string subFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"IndividualOutputTest");
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
            Assert.IsTrue(expectedFiles.Except(writtenFiles).Count() == 0);

            //check the second one is compressed and contains all the information
            writtenFiles = new HashSet<string>(Directory.GetFiles(Path.Combine(outputFolder, "compressFiles")).Select(v => Path.GetFileName(v)));
            //check the zip exists
            Assert.IsTrue(writtenFiles.Contains("Individual File Results.zip"));
            //check the original folder does not exist
            string[] subfolders = Directory.GetDirectories(Path.Combine(outputFolder, "compressFiles"));
            Assert.IsTrue(subfolders.Length == 0);
            ZipFile.ExtractToDirectory(Path.Combine(outputFolder, "compressFiles", "Individual File Results.zip"), Path.Combine(outputFolder, "compressFiles", "Individual File Results"));
            //read the extracted files
            writtenFiles = new HashSet<string>(Directory.GetFiles(Path.Combine(outputFolder, "compressFiles", "Individual File Results")).Select(v => Path.GetFileName(v).Substring(7)));
            //check they're the same
            Assert.IsTrue(expectedFiles.Except(writtenFiles).Count() == 0);

            //check the last one to make sure nothing was written except for the mzID files
            writtenFiles = new HashSet<string>(Directory.GetFiles(Path.Combine(outputFolder, "noFiles", "Individual File Results")).Select(v => Path.GetFileName(v).Substring(7)));
            Assert.IsTrue(writtenFiles.Count == 1);
            Assert.IsTrue(writtenFiles.Contains(".mzID"));

            Directory.Delete(outputFolder, true);

            //Do a check that we don't crash if there's only one file but somebody tries to zip the individual file results
            SearchTask weirdTask = new SearchTask();
            weirdTask.SearchParameters.CompressIndividualFiles = true;
            weirdTask.SearchParameters.WriteMzId = false;
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("weird", weirdTask) }, new List<string> { pathOne }, new List<DbForTask> { db }, outputFolder).Run();
            //check that a zip was not created
            writtenFiles = new HashSet<string>(Directory.GetFiles(Path.Combine(outputFolder, "weird")));
            Assert.IsFalse(writtenFiles.Contains("Individual File Results.zip"));
            Directory.Delete(subFolder, true);
        }

        [Test]
        public static void TestOpairFileOutput()
        {
            var task = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/GlycoSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);

            string currentDirecotry = Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData");

            Directory.CreateDirectory(currentDirecotry);
            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/P16150.fasta"), false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) }, new List<string> { spectraFile }, new List<DbForTask> { db }, currentDirecotry).Run();

            string outputDirectory = Path.Combine(currentDirecotry, @"Task");

            HashSet<string> writtenFiles = new HashSet<string>(Directory.GetFiles(outputDirectory).Select(v => Path.GetFileName(v)));

            //check that the first task wrote everything fine
            HashSet<string> expectedFiles = new HashSet<string>
            {
                "oglyco.psmtsv",
                "AutoGeneratedManuscriptProse.txt",
                "protein_oglyco_localization.tsv",
                "results.txt",
                "seen_oglyco_localization.tsv"
            };
            
            //check they're the same
            Assert.IsTrue(expectedFiles.Except(writtenFiles).Count() == 0);

            Directory.Delete(currentDirecotry, true);
        }
    }
}