using EngineLayer.DatabaseLoading;
using Nett;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Writing calibrated spectra as mgf, and choosing between the two mgf variants (MetaMorpheus#1673).
    /// </summary>
    [TestFixture]
    internal static class MgfOutputTests
    {
        /// <summary>
        /// OutputFormat is a plain enum, so Nett round-trips it as a string with no converter -- but only
        /// as long as nobody registers one for it. Every value, so a new one cannot be added untested.
        /// </summary>
        [Test]
        [TestCase(SpectraFileOutputFormat.MzML)]
        [TestCase(SpectraFileOutputFormat.Mgf)]
        [TestCase(SpectraFileOutputFormat.MgfMs2Only)]
        public static void OutputFormatSurvivesATomlRoundTrip(SpectraFileOutputFormat format)
        {
            var task = new CalibrationTask();
            task.CalibrationParameters.OutputFormat = format;

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "outputFormat.toml");

            try
            {
                Toml.WriteFile(task, path, MetaMorpheusTask.tomlConfig);
                Assert.That(File.ReadAllLines(path), Has.Member("OutputFormat = \"" + format + "\""));

                var reread = Toml.ReadFile<CalibrationTask>(path, MetaMorpheusTask.tomlConfig);
                Assert.That(reread.CalibrationParameters.OutputFormat, Is.EqualTo(format));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public static void OutputFormatDefaultsToMzMLSoNothingChangesUnlessAsked()
        {
            Assert.That(new CalibrationTask().CalibrationParameters.OutputFormat,
                Is.EqualTo(SpectraFileOutputFormat.MzML));
        }

        /// <summary>
        /// End to end: the extension, the block contents and the specification-conformance of each variant.
        /// The next task in the chain is handed this file by extension, so writing mgf under a .mzML name
        /// would hand it to the wrong reader.
        /// </summary>
        [Test]
        [NonParallelizable]
        [TestCase(SpectraFileOutputFormat.MzML, ".mzML")]
        [TestCase(SpectraFileOutputFormat.Mgf, ".mgf")]
        [TestCase(SpectraFileOutputFormat.MgfMs2Only, ".mgf")]
        public static void CalibrationWritesTheChosenFormat(SpectraFileOutputFormat format, string expectedExtension)
        {
            string unitTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "MgfOutput_" + format);
            string outputFolder = Path.Combine(unitTestFolder, "TaskOutput");
            Directory.CreateDirectory(outputFolder);

            try
            {
                string spectraPath = Path.Combine(unitTestFolder, "sample.mzML");
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "SmallCalibratible_Yeast.mzML"),
                    spectraPath, true);
                string database = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "smalldb.fasta");

                var task = new CalibrationTask();
                task.CalibrationParameters.OutputFormat = format;
                task.RunTask(outputFolder, new List<DbForTask> { new(database, false) },
                    new List<string> { spectraPath }, "test");

                string written = Path.Combine(outputFolder, "sample-calib" + expectedExtension);
                Assert.That(File.Exists(written), Is.True, "calibrated file should be " + expectedExtension);

                if (format == SpectraFileOutputFormat.MzML)
                {
                    return;
                }

                string[] lines = File.ReadAllLines(written);
                int blocks = lines.Count(l => l == "BEGIN IONS");
                int ms1Blocks = lines.Count(l => l == "MSLEVEL=1");
                int pepmasses = lines.Count(l => l.StartsWith("PEPMASS="));

                Assert.That(blocks, Is.GreaterThan(0));

                if (format == SpectraFileOutputFormat.MgfMs2Only)
                {
                    Assert.That(ms1Blocks, Is.Zero, "MgfMs2Only must not write MS1 scans");
                    Assert.That(pepmasses, Is.EqualTo(blocks),
                        "every block needs the mandatory PEPMASS, or MSToolkit-based readers exit on it");
                }
                else
                {
                    Assert.That(ms1Blocks, Is.GreaterThan(0), "Mgf must write MS1 scans so precursors survive");
                    Assert.That(pepmasses, Is.EqualTo(blocks - ms1Blocks));
                }
            }
            finally
            {
                Directory.Delete(unitTestFolder, true);
            }
        }
    }
}
