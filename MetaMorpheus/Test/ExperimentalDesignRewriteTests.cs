using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Reflection;
using EngineLayer;
using EngineLayer.Util;
using MassSpectrometry;
using NUnit.Framework;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Covers the experimental design that calibration and averaging rewrite for the next task. Every
    /// row of it has to name the file the task actually wrote, or quantification looks for something
    /// that is not there and skips - which is what #2192 reported.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class ExperimentalDesignRewriteTests
    {
        private static readonly List<string> Roots = new();

        [OneTimeTearDown]
        public static void OneTimeTearDown()
        {
            foreach (string root in Roots.Where(Directory.Exists))
            {
                Directory.Delete(root, true);
            }
        }

        /// <summary>
        /// WriteNewExperimentalDesignFile is private; reaching it directly is much cheaper than driving
        /// a calibration failure end to end, and it is the unit whose contract is at issue.
        /// </summary>
        private static void RewriteDesign(string oldDesignPath, string outputFolder,
            List<string> currentRawFileList, List<string> unsuccessful, Type task = null)
        {
            task ??= typeof(CalibrationTask);
            var method = task.GetMethod("WriteNewExperimentalDesignFile",
                BindingFlags.NonPublic | BindingFlags.Static);
            Assert.That(method, Is.Not.Null, task.Name + ".WriteNewExperimentalDesignFile was renamed");
            method.Invoke(null, new object[] { oldDesignPath, outputFolder, currentRawFileList, unsuccessful });
        }

        private static (string InDir, string OutDir) MakeDirs(string name)
        {
            string root = Path.Combine(TestContext.CurrentContext.TestDirectory, name);
            if (Directory.Exists(root)) Directory.Delete(root, true);
            Roots.Add(root);
            string inDir = Path.Combine(root, "in");
            string outDir = Path.Combine(root, "out");
            Directory.CreateDirectory(inDir);
            Directory.CreateDirectory(outDir);
            return (inDir, outDir);
        }

        /// <summary>The (FileName, Condition) rows of the rewritten design, in order.</summary>
        private static List<(string File, string Condition)> DesignRows(string outputFolder)
        {
            string path = Path.Combine(outputFolder, GlobalVariables.ExperimentalDesignFileName);
            Assert.That(File.Exists(path), "no experimental design was written into " + outputFolder);
            return File.ReadAllLines(path).Skip(1)
                .Where(l => !string.IsNullOrWhiteSpace(l))
                .Select(l => l.Split('\t'))
                .Select(cells => (cells[0], cells[1])).ToList();
        }

        private static List<string> DesignFileNames(string outputFolder) =>
            DesignRows(outputFolder).Select(r => r.File).ToList();

        [Test]
        public static void OneFileNameBeingAPrefixOfAnotherDoesNotPairThemTogether()
        {
            var (inDir, outDir) = MakeDirs("DesignRewritePrefix");

            // "run_1" is a substring of "run_10", which is the whole point.
            string calibratable = Path.Combine(inDir, "run_1.mzML");
            string failed = Path.Combine(inDir, "run_10.mzML");
            File.WriteAllText(calibratable, "x");
            File.WriteAllText(failed, "x");

            string oldDesign = ExperimentalDesign.WriteExperimentalDesignToFile(new List<SpectraFileInfo>
            {
                new SpectraFileInfo(calibratable, "CondA", 0, 0, 0),
                new SpectraFileInfo(failed, "CondB", 0, 0, 0)
            });

            // run_10 failed calibration, so its uncalibrated copy sits in the output folder.
            RewriteDesign(oldDesign, outDir, new List<string> { calibratable, failed },
                new List<string> { Path.Combine(outDir, "run_10.mzML") });

            var rows = DesignRows(outDir);
            Assert.That(rows, Is.EqualTo(new List<(string, string)>
            {
                ("run_1-calib.mzML", "CondA"),
                ("run_10.mzML", "CondB"),
            }), "each row must name its own file under its own condition: " + string.Join(", ", rows));
        }

        [Test]
        public static void ADesignRowNamesTheTrimmedPathTheRunLoopActuallyWrote()
        {
            var (inDir, outDir) = MakeDirs("DesignRewriteLongName");

            // Long enough that the calibrated path trips PathSafety's cap, as it does in the run loop.
            string longFile = Path.Combine(inDir, new string('a', 250) + ".mzML");
            File.WriteAllText(longFile, "x");

            string oldDesign = ExperimentalDesign.WriteExperimentalDesignToFile(
                new List<SpectraFileInfo> { new SpectraFileInfo(longFile, "CondA", 0, 0, 0) });

            // What CalibrationTask's run loop writes for this file.
            string writtenByRunLoop = Path.GetFileName(
                Path.Combine(outDir, Path.GetFileNameWithoutExtension(longFile) + "-calib" + ".mzML")
                    .ToSafeOutputPath("-calib.mzML"));
            Assert.That(writtenByRunLoop.Length, Is.LessThan(250),
                "premise: this name must actually be trimmed, or the test proves nothing");

            RewriteDesign(oldDesign, outDir, new List<string> { longFile }, new List<string>());

            Assert.That(DesignFileNames(outDir).Single(), Is.EqualTo(writtenByRunLoop));
        }

        /// <summary>
        /// GlobalVariables.GetFilenameWithoutExtension("exp.dia.d") is "expia": it removes every ".d".
        /// The skipped file is recorded in its input-path form, as CanFileTypeBeCalibrated does.
        /// </summary>
        [Test]
        public static void ASkippedDotDWhoseStemContainsDotDKeepsItsOwnRow()
        {
            var (inDir, outDir) = MakeDirs("DesignRewriteDotD");

            string calibratable = Path.Combine(inDir, "other.mzML");
            string skipped = Path.Combine(inDir, "exp.dia.d");

            string oldDesign = ExperimentalDesign.WriteExperimentalDesignToFile(new List<SpectraFileInfo>
            {
                new SpectraFileInfo(calibratable, "CondA", 0, 0, 0),
                new SpectraFileInfo(skipped, "CondB", 0, 0, 0)
            });

            RewriteDesign(oldDesign, outDir, new List<string> { calibratable, skipped }, new List<string> { skipped });

            Assert.That(DesignRows(outDir), Is.EqualTo(new List<(string, string)>
            {
                ("other-calib.mzML", "CondA"),
                ("exp.dia.d", "CondB"),
            }));
        }

        [Test]
        public static void ASkippedFileDoesNotClaimTheRowOfACalibratedFileWithTheSameStem()
        {
            var (inDir, outDir) = MakeDirs("DesignRewriteSameStem");

            string calibratable = Path.Combine(inDir, "run.mzML");
            string skipped = Path.Combine(inDir, "run.mgf");

            string oldDesign = ExperimentalDesign.WriteExperimentalDesignToFile(new List<SpectraFileInfo>
            {
                new SpectraFileInfo(calibratable, "CondA", 0, 0, 0),
                new SpectraFileInfo(skipped, "CondB", 0, 0, 0)
            });

            RewriteDesign(oldDesign, outDir, new List<string> { calibratable, skipped }, new List<string> { skipped });

            Assert.That(DesignRows(outDir), Is.EqualTo(new List<(string, string)>
            {
                ("run-calib.mzML", "CondA"),
                ("run.mgf", "CondB"),
            }));
        }

        /// <summary>
        /// A skipped file's row keeps its input path. Listed first, it used to send the rewritten design
        /// into the input folder, over the user's own design.
        /// </summary>
        [Test]
        [TestCase(typeof(CalibrationTask), "-calib")]
        [TestCase(typeof(SpectralAveragingTask), "-averaged")]
        public static void ASkippedFileListedFirstDoesNotOverwriteTheInputDesign(Type task, string suffix)
        {
            var (inDir, outDir) = MakeDirs("DesignRewriteSkippedFirst" + task.Name);

            string skipped = Path.Combine(inDir, "sample2.mgf");
            string processed = Path.Combine(inDir, "sample1.mzML");

            string oldDesign = ExperimentalDesign.WriteExperimentalDesignToFile(new List<SpectraFileInfo>
            {
                new SpectraFileInfo(skipped, "CondB", 0, 0, 0),
                new SpectraFileInfo(processed, "CondA", 0, 0, 0)
            });
            string inputDesignBefore = File.ReadAllText(oldDesign);

            RewriteDesign(oldDesign, outDir, new List<string> { skipped, processed }, new List<string> { skipped }, task);

            Assert.That(File.ReadAllText(oldDesign), Is.EqualTo(inputDesignBefore), "the input design was overwritten");
            Assert.That(DesignRows(outDir), Is.EqualTo(new List<(string, string)>
            {
                ("sample2.mgf", "CondB"),
                ("sample1" + suffix + ".mzML", "CondA"),
            }));
        }

        /// <summary>
        /// GlobalVariables.GetFilenameWithoutExtension("exp.mzML.v2.mzML") is "exp.v2"; the run loop
        /// writes exp.mzML.v2-averaged.mzML.
        /// </summary>
        [Test]
        public static void AnAveragedRowNamesTheFileTheRunLoopWroteWhenTheExtensionRecursInTheStem()
        {
            var (inDir, outDir) = MakeDirs("DesignRewriteAveragedStem");

            string file = Path.Combine(inDir, "exp.mzML.v2.mzML");

            string oldDesign = ExperimentalDesign.WriteExperimentalDesignToFile(
                new List<SpectraFileInfo> { new SpectraFileInfo(file, "CondA", 0, 0, 0) });

            RewriteDesign(oldDesign, outDir, new List<string> { file }, new List<string>(), typeof(SpectralAveragingTask));

            Assert.That(DesignRows(outDir), Is.EqualTo(new List<(string, string)>
            {
                ("exp.mzML.v2-averaged.mzML", "CondA"),
            }));
        }
    }
}
