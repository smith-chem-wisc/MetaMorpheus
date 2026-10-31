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
    /// Covers the experimental design that calibration rewrites for the next task. Every row of it has
    /// to name the file calibration actually wrote, or quantification looks for something that is not
    /// there and skips - which is what #2192 reported.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class ExperimentalDesignRewriteTests
    {
        /// <summary>
        /// WriteNewExperimentalDesignFile is private; reaching it directly is much cheaper than driving
        /// a calibration failure end to end, and it is the unit whose contract is at issue.
        /// </summary>
        private static void RewriteDesign(string oldDesignPath, string outputFolder,
            List<string> currentRawFileList, List<string> unsuccessfullyCalibrated)
        {
            var method = typeof(CalibrationTask).GetMethod("WriteNewExperimentalDesignFile",
                BindingFlags.NonPublic | BindingFlags.Static);
            Assert.That(method, Is.Not.Null, "CalibrationTask.WriteNewExperimentalDesignFile was renamed");
            method.Invoke(null, new object[] { oldDesignPath, outputFolder, currentRawFileList, unsuccessfullyCalibrated });
        }

        private static (string InDir, string OutDir) MakeDirs(string name)
        {
            string root = Path.Combine(TestContext.CurrentContext.TestDirectory, name);
            if (Directory.Exists(root)) Directory.Delete(root, true);
            string inDir = Path.Combine(root, "in");
            string outDir = Path.Combine(root, "out");
            Directory.CreateDirectory(inDir);
            Directory.CreateDirectory(outDir);
            return (inDir, outDir);
        }

        /// <summary>The FileName column of the rewritten design, in order.</summary>
        private static List<string> DesignFileNames(string outputFolder)
        {
            string path = Path.Combine(outputFolder, GlobalVariables.ExperimentalDesignFileName);
            Assert.That(File.Exists(path), "no experimental design was written");
            return File.ReadAllLines(path).Skip(1)
                .Where(l => !string.IsNullOrWhiteSpace(l))
                .Select(l => l.Split('\t')[0]).ToList();
        }

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

            var names = DesignFileNames(outDir);

            Assert.That(names, Has.Count.EqualTo(2));
            Assert.That(names, Is.Unique, "a file was named twice, so one row lost its own data: " + string.Join(", ", names));
            Assert.That(names, Does.Contain("run_1-calib.mzML"),
                "run_1 calibrated successfully, so its row must name its calibrated file: " + string.Join(", ", names));
            Assert.That(names, Does.Contain("run_10.mzML"));
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
    }
}
