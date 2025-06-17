using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class MsDataFileTest
    {
        [OneTimeSetUp]
        public static void Setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;
        }

        [Test]
        [TestCase(@"TestData\ok.mgf", @"TestData\okk.xml")]
        [TestCase(@"TestData\snippet.d", @"TestData\gapdh.fasta")]
        public static void TestQuantificationDoesntCrashOnUnsupportedFiles(string filepath, string dbPath)
        {
            //The purpose of this test is to ensure that mgfs and timsTOF (.d) files can be run without crashing.

            //Whenever a new feature is added that may require things an mgf does not have,
            //there should be a check that prevents mgfs from using that feature.
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestLoadAndRunMgf");

            SearchTask task1 = new()
            {
                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    DoLabelFreeQuantification = true
                }
            };
            List<(string, MetaMorpheusTask)> taskList = new()
            {
                ("task1", task1),
            };
            //run!

            var engine = new EverythingRunnerEngine(taskList, new List<string> { filepath }, new List<DbForTask> { new DbForTask(dbPath, false) }, outputFolder);
            engine.Run();
            //Just don't crash! There should also be at least one psm at 1% FDR, but can't check for that.
            Directory.Delete(outputFolder, true);
        }

        [Test]
        [TestCase(@"TestData\ok.mgf", @"TestData\okk.xml")]
        [TestCase(@"TestData\snippet.d", @"TestData\gapdh.fasta")]
        public static void TestCalibrationDoesntCrashOnUnsupportedFiles(string filepath, string dbPath)
        {
            //The purpose of this test is to ensure that mgfs and timsTOF (.d) files can be run without crashing.

            //Whenever a new feature is added that may require things an mgf does not have,
            //there should be a check that prevents mgfs from using that feature.
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestCalibrateAndSearchMgf");

            CalibrationTask task1 = new();
            SearchTask task2 = new()
            {
                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    DoLabelFreeQuantification = false
                }
            };
            List<(string, MetaMorpheusTask)> taskList = new()
            {
                ("task1", task1),
                ("task2", task2)
            };
            //run!

            var engine = new EverythingRunnerEngine(taskList, new List<string> { filepath }, new List<DbForTask> { new DbForTask(dbPath, false) }, outputFolder);
            engine.Run();
            //Just don't crash! There should also be at least one psm at 1% FDR, but can't check for that.
            Directory.Delete(outputFolder, true);
        }

        [Test]
        [TestCase(@"TestData\ok.mgf", @"TestData\okk.xml")]
        [TestCase(@"TestData\snippet.d", @"TestData\gapdh.fasta")]
        public static void TestAveragingDoesntCrashOnUnsupportedFiles(string filepath, string dbPath)
        {
            //The purpose of this test is to ensure that mgfs and timsTOF (.d) files can be run without crashing.

            //Whenever a new feature is added that may require things an mgf does not have,
            //there should be a check that prevents mgfs from using that feature.
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestCalibrateAndSearchMgf");

            SpectralAveragingTask task1 = new()
            {
                CommonParameters = new CommonParameters()
            };
            SearchTask task2 = new()
            {
                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    DoLabelFreeQuantification = false
                }
            };
            List<(string, MetaMorpheusTask)> taskList = new()
            {
                ("task1", task1),
                ("task2", task2)
            };
            //run!

            var engine = new EverythingRunnerEngine(taskList, new List<string> { filepath }, new List<DbForTask> { new DbForTask(dbPath, false) }, outputFolder);
            engine.Run();
            //Just don't crash! There should also be at least one psm at 1% FDR, but can't check for that.
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestCompressionDecompression()
        {
            string testInputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"CompressionTest");
            DirectoryInfo testDirectory = new(testInputFolder);
            MyFileManager.CompressDirectory(testDirectory);

            foreach (FileInfo file in testDirectory.GetFiles())
            {
                Assert.That(file.Extension, Is.EqualTo(".gz"));
            }

            MyFileManager.DecompressDirectory(testDirectory);

            foreach (FileInfo file in testDirectory.GetFiles())
            {
                Assert.That(file.Extension, Is.Not.EqualTo(".gz"));
            }
        }

        [Test]
        public static void TestMs2ScanWithSpecificMass()
        {
            Ms2ScanWithSpecificMass scanB = new(
                new MsDataScan(
                    new MzSpectrum(Array.Empty<double>(), Array.Empty<double>(), false),
                    2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null),
                100, 1, null, new CommonParameters(), null);

            var closestExperimentalMassB = scanB.GetClosestExperimentalIsotopicEnvelope(10);

            Assert.That(closestExperimentalMassB, Is.Null);
        }
    }
}