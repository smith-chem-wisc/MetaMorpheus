using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.Util;
using NUnit.Framework;
using SpectralAveraging;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// EverythingRunnerEngine hands each task the files the task before it wrote, plus every input that task wrote nothing
    /// for: one it skips (.mgf, .d, .msalign) or could not calibrate. It found those inputs by stripping both -calib and
    /// -averaged from each written name and looking for what was left among the inputs. An input that was itself a -calib
    /// file was never found that way (x-calib-averaged strips to x, not x-calib), and neither was a name trimmed to the
    /// path cap, so the input was handed on beside the file written for it and the search read it twice (#2942). When no
    /// input is skipped every input is replaced, so each test here has a skipped .mgf.
    /// </summary>
    [TestFixture]
    public static class SpectraHandOffTests
    {
        private static string TestData(string name) => Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", name);

        private static List<DbForTask> YeastDatabase() => new() { new DbForTask(TestData("smalldb.fasta"), false) };

        private static string Fresh(string name)
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, name);
            if (Directory.Exists(folder)) Directory.Delete(folder, true);
            Directory.CreateDirectory(Path.Combine(folder, "TaskOutput"));
            return folder;
        }

        private static MetaMorpheusTask DerivingTask(string kind) => kind == "calibrate"
            ? new CalibrationTask()
            : new SpectralAveragingTask(new SpectralAveragingParameters { SpectraFileAveragingType = SpectraFileAveragingType.AverageDdaScans, NumberOfScansToAverage = 5 })
            {
                CommonParameters = new CommonParameters(taskDescriptor: "AveragingTask", maxThreadsToUsePerFile: 2)
            };

        /// <summary>The files the search was given, in the order it was given them.</summary>
        private static List<string> SearchedFiles(SearchTask search) =>
            search.FileSpecificParameters.Select(file => Path.GetFullPath(file.FileName)).ToList();

        /// <summary>
        /// The reproduction in #2942: Calibrate -> Average -> Search with an .mgf that both tasks skip searched the
        /// -calib-averaged file, the -calib file it was averaged from, and the .mgf. Calibrate -> Calibrate takes the same
        /// path. The search must get the last task's file and the .mgf, once each.
        /// </summary>
        [TestCase("calibrate", "average", "filename1-calib-averaged.mzML")]
        [TestCase("calibrate", "calibrate", "filename1-calib-calib.mzML")]
        public static void SearchAfterTwoTasks_ReadsEachFileOnce_WhenAnInputIsSkipped(string first, string second, string lastWritten)
        {
            string folder = Fresh("SpectraHandOffTwoTasks");
            string outputFolder = Path.Combine(folder, "TaskOutput");
            string spectraFile = Path.Combine(folder, "filename1.mzML");
            File.Copy(TestData("SmallCalibratible_Yeast.mzML"), spectraFile, true);
            string mgf = Path.Combine(folder, "tinySemi.mgf");
            File.Copy(TestData("tinySemi.mgf"), mgf, true);

            var search = new SearchTask();
            var tasks = new List<(string, MetaMorpheusTask)> { ("first", DerivingTask(first)), ("second", DerivingTask(second)), ("search", search) };
            try
            {
                new EverythingRunnerEngine(tasks, new List<string> { spectraFile, mgf }, YeastDatabase(), outputFolder).Run();

                string written = Path.Combine(outputFolder, "second", lastWritten);
                Assert.That(File.Exists(written), "the fixture only means anything if both tasks write a file for the mzML");
                Assert.That(SearchedFiles(search), Is.EqualTo(new[] { Path.GetFullPath(written), Path.GetFullPath(mgf) }));
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// The run loop names a -calib file through ToSafeOutputPath, so a -calib path over the path cap is written trimmed.
        /// Stripping -calib from the trimmed name did not give the input's name, so Calibrate -> Search with a skipped .mgf
        /// searched the trimmed -calib file and the uncalibrated input both. The input path stays under the cap; only the
        /// -calib path is over it.
        /// </summary>
        [Test]
        public static void SearchAfterCalibration_ReadsEachFileOnce_WhenTheCalibratedNameIsTrimmed()
        {
            string folder = Fresh("SpectraHandOffTrimmed");
            string outputFolder = Path.Combine(folder, "TaskOutput");
            string calibrationFolder = Path.Combine(outputFolder, "calibrate");
            string calibratedEnding = CalibrationTask.CalibSuffix + ".mzML";
            // five characters over the cap once the -calib ending is added in the calibration folder
            int nameLength = PathSafety.DefaultMaxPath - (calibrationFolder.Length + 1) - calibratedEnding.Length + 5;
            string spectraFile = Path.Combine(folder, new string('a', nameLength) + ".mzML");
            Assert.That(spectraFile.Length, Is.LessThan(PathSafety.DefaultMaxPath),
                "the fixture only means anything if the input path itself is under the cap");
            File.Copy(TestData("SmallCalibratible_Yeast.mzML"), spectraFile, true);
            string mgf = Path.Combine(folder, "tinySemi.mgf");
            File.Copy(TestData("tinySemi.mgf"), mgf, true);

            // per-file outputs would be named for the untrimmed name, past the cap
            var search = new SearchTask();
            search.SearchParameters.WriteIndividualFiles = false;
            var tasks = new List<(string, MetaMorpheusTask)> { ("calibrate", new CalibrationTask()), ("search", search) };
            try
            {
                new EverythingRunnerEngine(tasks, new List<string> { spectraFile, mgf }, YeastDatabase(), outputFolder).Run();

                var calibrated = Directory.GetFiles(calibrationFolder, "*" + calibratedEnding);
                Assert.That(calibrated, Has.Length.EqualTo(1), "the fixture only means anything if calibration succeeds");
                Assert.That(Path.GetFileName(calibrated[0]).Length, Is.LessThan(nameLength + calibratedEnding.Length),
                    "and if the run loop trimmed the calibrated name");
                Assert.That(SearchedFiles(search), Is.EqualTo(new[] { Path.GetFullPath(calibrated[0]), Path.GetFullPath(mgf) }));
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// The TMT half of #2942. Calibration re-writes TmtDesign.txt for the files it wrote, and the search reads the design
        /// beside its first file and refuses one that does not name every file it searches. Handed the -calib file a second
        /// time, the search after Calibrate -> Calibrate warned "The TMT design did not contain the file(s)" and skipped
        /// multiplex quantification. The design and the search's files must agree. (Averaging re-writes TmtDesign.txt only
        /// from #2938 on, so Calibrate -> Average stands in for itself only once that merges.) The yeast file carries no
        /// reporter ions; the design is what is tested, and the search archives it into its output only after reading it.
        /// </summary>
        [Test]
        public static void SearchAfterCalibratingTwice_ReadsTheTmtDesign_WhenAnInputIsSkipped()
        {
            string folder = Fresh("SpectraHandOffTmt");
            string outputFolder = Path.Combine(folder, "TaskOutput");
            string spectraFile = Path.Combine(folder, "filename1.mzML");
            File.Copy(TestData("SmallCalibratible_Yeast.mzML"), spectraFile, true);
            string mgf = Path.Combine(folder, "tinySemi.mgf");
            File.Copy(TestData("tinySemi.mgf"), mgf, true);
            CalibrationTests.WriteTmtDesignBeside(spectraFile, mgf);

            var search = new SearchTask();
            search.SearchParameters.DoMultiplexQuantification = true;
            search.SearchParameters.MultiplexModId = "TMT11";
            search.SearchParameters.DoParsimony = true;
            var tasks = new List<(string, MetaMorpheusTask)> { ("first", new CalibrationTask()), ("second", new CalibrationTask()), ("search", search) };

            var warnings = new List<string>();
            EventHandler<StringEventArgs> onWarn = (_, e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += onWarn;
            try
            {
                new EverythingRunnerEngine(tasks, new List<string> { spectraFile, mgf }, YeastDatabase(), outputFolder).Run();

                Assert.That(File.Exists(Path.Combine(outputFolder, "second", "filename1-calib-calib.mzML")),
                    "the fixture only means anything if both calibrations succeed");
                List<string> searched = SearchedFiles(search);
                var design = TmtExperimentalDesign.Read(Path.Combine(Path.GetDirectoryName(searched[0]), GlobalVariables.TmtExperimentalDesignFileName),
                    searched, out var errors);
                Assert.Multiple(() =>
                {
                    Assert.That(errors, Is.Empty);
                    Assert.That(design.Select(file => Path.GetFullPath(file.FullFilePathWithExtension).ToLowerInvariant()),
                        Is.EquivalentTo(searched.Select(file => file.ToLowerInvariant())));
                    Assert.That(warnings, Has.None.Contains("TMT design"));
                    Assert.That(File.Exists(Path.Combine(outputFolder, "search", GlobalVariables.TmtExperimentalDesignFileName)), Is.True,
                        string.Join(Environment.NewLine, warnings));
                });
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= onWarn;
                Directory.Delete(folder, true);
            }
        }
    }
}
