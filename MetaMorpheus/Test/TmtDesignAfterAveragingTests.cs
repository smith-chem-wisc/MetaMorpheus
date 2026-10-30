using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using NUnit.Framework;
using SpectralAveraging;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Averaging re-wrote ExperimentalDesign.tsv for its -averaged files but not TmtDesign.txt. The search looks for the
    /// TMT design beside the files it searches, so a TMT search after averaging warned "No TmtDesign.txt found next to the
    /// spectra files" and wrote no quantification table. Averaging now re-writes it as calibration does (#2904), through
    /// the same <see cref="MetaMorpheusTask"/> method.
    /// </summary>
    [TestFixture]
    public static class TmtDesignAfterAveragingTests
    {
        private static string Fresh(string name)
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, name);
            if (Directory.Exists(folder)) Directory.Delete(folder, true);
            Directory.CreateDirectory(Path.Combine(folder, "TaskOutput"));
            return folder;
        }

        private static string CopyOfYeast(string folder, string name)
        {
            string spectraFile = Path.Combine(folder, name);
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "SmallCalibratible_Yeast.mzML"), spectraFile, true);
            return spectraFile;
        }

        private static SpectralAveragingTask AveragingTask() =>
            new(new SpectralAveragingParameters { SpectraFileAveragingType = SpectraFileAveragingType.AverageDdaScans, NumberOfScansToAverage = 5 })
            {
                CommonParameters = new CommonParameters(taskDescriptor: "AveragingTask", maxThreadsToUsePerFile: 2)
            };

        private static List<DbForTask> YeastDatabase() =>
            new() { new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "smalldb.fasta"), false) };

        [Test]
        public static void TmtDesignIsRewrittenForTheAveragedFile()
        {
            string folder = Fresh("TmtDesignAveragingTest");
            string outputFolder = Path.Combine(folder, "TaskOutput");
            string spectraFile = CopyOfYeast(folder, "filename1.mzML");
            CalibrationTests.WriteTmtDesignBeside(spectraFile);
            string inputDesign = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);
            var original = TmtExperimentalDesign.Read(inputDesign, new List<string> { spectraFile }, out var originalErrors);
            Assert.That(originalErrors, Is.Empty);
            byte[] inputBytes = File.ReadAllBytes(inputDesign);

            try
            {
                AveragingTask().RunTask(outputFolder, new List<DbForTask>(), new List<string> { spectraFile }, "test");

                string averagedFile = Path.Combine(outputFolder, "filename1-averaged.mzML");
                Assert.That(File.Exists(averagedFile), "the fixture only means anything if averaging succeeds");

                var rewritten = TmtExperimentalDesign.Read(Path.Combine(outputFolder, GlobalVariables.TmtExperimentalDesignFileName),
                    new List<string> { averagedFile }, out var errors);

                Assert.That(errors, Is.Empty);
                Assert.Multiple(() =>
                {
                    Assert.That(Path.GetFullPath(rewritten.Single().FullFilePathWithExtension), Is.EqualTo(Path.GetFullPath(averagedFile)).IgnoreCase);
                    Assert.That(rewritten[0].Plex, Is.EqualTo(original[0].Plex));
                    Assert.That(rewritten[0].Fraction, Is.EqualTo(original[0].Fraction));
                    Assert.That(rewritten[0].TechnicalReplicate, Is.EqualTo(original[0].TechnicalReplicate));
                    Assert.That(rewritten[0].Annotations.Select(a => (a.Tag, a.SampleName, a.Condition, a.BiologicalReplicate, a.SampleType)),
                        Is.EqualTo(original[0].Annotations.Select(a => (a.Tag, a.SampleName, a.Condition, a.BiologicalReplicate, a.SampleType))));
                    Assert.That(File.ReadAllBytes(inputDesign), Is.EqualTo(inputBytes), "the input design is the user's");
                });
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// Averaging skips an .mgf and writes nothing for it, so the task after it searches the original
        /// (EverythingRunnerEngine passes on every file missing from NewSpectra), and the re-written design names it so.
        /// </summary>
        [Test]
        public static void TmtDesignNamesTheOriginalOfAFileAveragingSkips()
        {
            string folder = Fresh("TmtDesignAveragingSkipTest");
            string outputFolder = Path.Combine(folder, "TaskOutput");
            string spectraFile = CopyOfYeast(folder, "filename1.mzML");
            string mgf = Path.Combine(folder, "skipped.mgf");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "ok.mgf"), mgf, true);
            CalibrationTests.WriteTmtDesignBeside(spectraFile, mgf);

            try
            {
                AveragingTask().RunTask(outputFolder, new List<DbForTask>(), new List<string> { spectraFile, mgf }, "test");

                string averagedFile = Path.Combine(outputFolder, "filename1-averaged.mzML");
                Assert.That(File.Exists(averagedFile), "the fixture only means anything if the mzML is averaged");
                Assert.That(File.Exists(Path.Combine(outputFolder, "skipped-averaged.mzML")), Is.False);

                var rewritten = TmtExperimentalDesign.Read(Path.Combine(outputFolder, GlobalVariables.TmtExperimentalDesignFileName),
                    new List<string> { averagedFile, mgf }, out var errors);

                Assert.That(errors, Is.Empty);
                Assert.That(rewritten.Select(f => Path.GetFullPath(f.FullFilePathWithExtension).ToLowerInvariant()),
                    Is.EquivalentTo(new[] { Path.GetFullPath(averagedFile).ToLowerInvariant(), Path.GetFullPath(mgf).ToLowerInvariant() }));
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// End to end, the case measured on 2026-10-20: Average -> Search, and Calibrate -> Average -> Search, warned "No
        /// TmtDesign.txt found next to the spectra files" and wrote no quantification table. The yeast file carries no
        /// reporter ions; the design is what is tested, and the search archives it into its output only after reading it.
        /// </summary>
        [TestCase(false)]
        [TestCase(true)]
        public static void SearchAfterAveraging_FindsTheTmtDesign(bool calibrateFirst)
        {
            string folder = Fresh("TmtAverageThenSearchTest");
            string outputFolder = Path.Combine(folder, "TaskOutput");
            string spectraFile = CopyOfYeast(folder, "filename1.mzML");
            CalibrationTests.WriteTmtDesignBeside(spectraFile);

            var searchTask = new SearchTask();
            searchTask.SearchParameters.DoMultiplexQuantification = true;
            searchTask.SearchParameters.MultiplexModId = "TMT11";
            searchTask.SearchParameters.DoParsimony = true;
            var tasks = new List<(string, MetaMorpheusTask)>();
            if (calibrateFirst) tasks.Add(("calibrate", new CalibrationTask()));
            tasks.Add(("average", AveragingTask()));
            tasks.Add(("search", searchTask));

            var warnings = new List<string>();
            EventHandler<StringEventArgs> onWarn = (_, e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += onWarn;
            try
            {
                new EverythingRunnerEngine(tasks, new List<string> { spectraFile }, YeastDatabase(), outputFolder).Run();

                string averaged = calibrateFirst ? "filename1-calib-averaged.mzML" : "filename1-averaged.mzML";
                Assert.That(File.Exists(Path.Combine(outputFolder, "average", averaged)),
                    "the fixture only means anything if every step before the search writes its file");
                Assert.That(warnings, Has.None.Contains("No " + GlobalVariables.TmtExperimentalDesignFileName + " found"));
                Assert.That(File.Exists(Path.Combine(outputFolder, "search", GlobalVariables.TmtExperimentalDesignFileName)), Is.True,
                    string.Join(Environment.NewLine, warnings));
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= onWarn;
                Directory.Delete(folder, true);
            }
        }
    }
}
