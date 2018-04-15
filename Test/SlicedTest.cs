using Nett;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class SlicedTest
    {
        #region Public Methods

        [Test]
        public static void SlicedTest1()
        {
            var task = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"SlicedSearchTaskConfig.toml"), MetaMorpheusTask.tomlConfig);

            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-db.fasta"), false);
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-raw.mzML");
            EverythingRunnerEngine a = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) }, new List<string> { raw }, new List<DbForTask> { db }, Environment.CurrentDirectory);

            a.Run();

            var thisTaskOutputFolder = MySetUpClass.outputFolder;

            var peaks = Path.Combine(thisTaskOutputFolder, "Task", "sliced-raw_QuantifiedPeaks.tsv");

            Assert.AreEqual(2, File.ReadLines(peaks).Count());

            var psms = Path.Combine(thisTaskOutputFolder, "Task", "sliced-raw_PSMs.psmtsv");

            Assert.AreEqual(3, File.ReadLines(psms).Count());
            var protGroups = Path.Combine(thisTaskOutputFolder, "Task", "sliced-raw_ProteinGroups.tsv");

            Assert.AreEqual(2, File.ReadLines(protGroups).Count());
        }


        [Test]
        public static void FaFormatTest()
        {
            var task = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"SlicedSearchTaskConfig.toml"), MetaMorpheusTask.tomlConfig);

            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-db.fa"), false);
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-raw.mzML");
            EverythingRunnerEngine a = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) }, new List<string> { raw }, new List<DbForTask> { db }, Environment.CurrentDirectory);

            a.Run();

            var thisTaskOutputFolder = MySetUpClass.outputFolder;

            var peaks = Path.Combine(thisTaskOutputFolder, "Task", "sliced-raw_QuantifiedPeaks.tsv");

            Assert.AreEqual(2, File.ReadLines(peaks).Count());

            var psms = Path.Combine(thisTaskOutputFolder, "Task", "sliced-raw_PSMs.psmtsv");

            Assert.AreEqual(3, File.ReadLines(psms).Count());
            var protGroups = Path.Combine(thisTaskOutputFolder, "Task", "sliced-raw_ProteinGroups.tsv");

            Assert.AreEqual(2, File.ReadLines(protGroups).Count());
        }

        [Test]
        public static void CalibrationTest()
        {
            var preCalibrationSearch = new SearchTask();
            var calibration = new CalibrationTask();
            var postCalibrationSearch = new SearchTask();
            var gptmd = new GptmdTask();
            var postGptmdSearch = new SearchTask();

            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallYeastDb.fasta"), false);
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratibleYeast.mzml");

            EverythingRunnerEngine a = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)>
            { ("Task1", preCalibrationSearch),
                ("Task2", calibration),
                ("Task3", postCalibrationSearch),
                ("Task4", gptmd),
                ("Task5", postGptmdSearch)
            }, 
                new List<string> { raw }, 
                new List<DbForTask> { db }, 
                Environment.CurrentDirectory);

            a.Run();

            //var thisTaskOutputFolder = MySetUpClass.outputFolder;

            //var peaks = Path.Combine(thisTaskOutputFolder, "Task", "sliced-raw_QuantifiedPeaks.tsv");

            //Assert.AreEqual(2, File.ReadLines(peaks).Count());

            //var psms = Path.Combine(thisTaskOutputFolder, "Task", "sliced-raw_PSMs.psmtsv");

            //Assert.AreEqual(3, File.ReadLines(psms).Count());
            //var protGroups = Path.Combine(thisTaskOutputFolder, "Task", "sliced-raw_ProteinGroups.tsv");

            //Assert.AreEqual(2, File.ReadLines(protGroups).Count());
        }

        #endregion Public Methods
    }
}