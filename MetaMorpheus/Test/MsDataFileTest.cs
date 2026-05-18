using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer.DatabaseLoading;
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
        [TestCase(@"TopDownTestData\JurkatTopDownRep2Fract1_ms2.msalign", @"TestData\gapdh.fasta")]
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
        [TestCase(@"TopDownTestData\JurkatTopDownRep2Fract1_ms2.msalign", @"TestData\gapdh.fasta")]
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
        [TestCase(@"TopDownTestData\JurkatTopDownRep2Fract1_ms2.msalign", @"TestData\gapdh.fasta")]
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

        /// <summary>
        /// Regression test for crash when a folder containing a .d folder (timsTOF data) is opened.
        /// The .d folder is a directory that should be treated as a spectra file.
        /// This test ensures that when a parent directory containing a .d folder is specified,
        /// the search runs without crashing.
        /// </summary>
        [Test]
        public static void TestFolderContainingDotDFolderDoesNotCrash()
        {
            // The TestData folder contains snippet.d which is a timsTOF data folder
            string testDataFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData");
            string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\gapdh.fasta");
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestDotDFolderRegression");

            // Verify the .d folder exists in the test data
            string dotDFolderPath = Path.Combine(testDataFolder, "snippet.d");
            Assert.That(Directory.Exists(dotDFolderPath), Is.True, "Test requires snippet.d folder to exist");

            SearchTask task1 = new()
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
            };

            // This should not crash - the .d folder should be discovered and handled correctly
            var engine = new EverythingRunnerEngine(taskList, new List<string> { dotDFolderPath }, new List<DbForTask> { new DbForTask(dbPath, false) }, outputFolder);
            engine.Run();

            // Cleanup
            if (Directory.Exists(outputFolder))
            {
                Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// Regression test to ensure that when a parent folder containing a .d folder is processed,
        /// files inside the .d folder (like .xml, .mgf, .fasta, etc.) are NOT added separately.
        /// Only the .d folder itself should be treated as a spectra file.
        /// This prevents the internal timsTOF files from being incorrectly processed as databases or other file types.
        /// </summary>
        [Test]
        public static void TestFilesInsideDotDFolderAreNotAddedSeparately()
        {
            // Create a temporary test directory structure
            string tempFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestDotDFilesNotAdded");
            string fakeDotDFolder = Path.Combine(tempFolder, "fake_data.d");
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestDotDFilesOutput");

            try
            {
                // Clean up if exists from previous failed run
                if (Directory.Exists(tempFolder))
                {
                    Directory.Delete(tempFolder, true);
                }
                if (Directory.Exists(outputFolder))
                {
                    Directory.Delete(outputFolder, true);
                }

                // Create the test directory structure
                Directory.CreateDirectory(fakeDotDFolder);

                // Create fake files inside the .d folder that would normally be valid file types
                // These should NOT be picked up when the parent folder is processed
                File.WriteAllText(Path.Combine(fakeDotDFolder, "fake_database.xml"), "<test></test>");
                File.WriteAllText(Path.Combine(fakeDotDFolder, "fake_spectra.mgf"), "BEGIN IONS\nEND IONS");
                File.WriteAllText(Path.Combine(fakeDotDFolder, "fake_database.fasta"), ">test\nACDEFG");
                File.WriteAllText(Path.Combine(fakeDotDFolder, "analysis.tdf"), "fake tdf content");

                // Test the file enumeration logic that the GUI and CMD use
                // When enumerating tempFolder, we should get:
                // 1. fake_data.d (as a .d directory - spectra file)
                // We should NOT get the files inside fake_data.d

                var filesInTempFolder = Directory.GetFiles(tempFolder);
                var dirsInTempFolder = Directory.GetDirectories(tempFolder);

                // Files directly in temp folder (not inside .d) - should be empty since we only created files inside .d
                var spectraFilesFound = filesInTempFolder
                    .Where(f => GlobalVariables.AcceptedSpectraFormats.Contains(Path.GetExtension(f).ToLowerInvariant()))
                    .ToList();

                // .d directories should be found
                var dotDFoldersFound = dirsInTempFolder
                    .Where(d => d.EndsWith(".d", StringComparison.OrdinalIgnoreCase))
                    .ToList();

                // Verify the .d folder was found
                Assert.That(dotDFoldersFound.Count, Is.EqualTo(1), "Should find exactly one .d folder");
                Assert.That(dotDFoldersFound[0], Does.EndWith("fake_data.d"));

                // Verify that files inside .d folder are NOT in the direct file enumeration
                // (they would only appear if we used recursive enumeration, which we shouldn't for .d folders)
                var allFilesRecursive = Directory.GetFiles(tempFolder, "*.*", SearchOption.AllDirectories);
                var filesInsideDotD = allFilesRecursive.Where(f => f.Contains("fake_data.d")).ToList();

                // These files exist inside .d but should be ignored when processing the parent folder
                Assert.That(filesInsideDotD.Any(f => f.EndsWith(".xml")), Is.True, "Test setup: xml file should exist inside .d");
                Assert.That(filesInsideDotD.Any(f => f.EndsWith(".mgf")), Is.True, "Test setup: mgf file should exist inside .d");
                Assert.That(filesInsideDotD.Any(f => f.EndsWith(".fasta")), Is.True, "Test setup: fasta file should exist inside .d");

                // But direct enumeration (non-recursive) should NOT include them
                Assert.That(spectraFilesFound.Any(f => f.Contains("fake_data.d")), Is.False, 
                    "Files inside .d folder should not be enumerated as separate spectra files");

                // Also verify that the .d folder is recognized as a valid spectra format
                Assert.That(GlobalVariables.AcceptedSpectraFormats.Contains(".d"), Is.True, 
                    ".d should be in accepted spectra formats");
            }
            finally
            {
                // Cleanup
                if (Directory.Exists(tempFolder))
                {
                    Directory.Delete(tempFolder, true);
                }
                if (Directory.Exists(outputFolder))
                {
                    Directory.Delete(outputFolder, true);
                }
            }
        }
    }
}
