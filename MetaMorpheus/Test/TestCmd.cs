using NUnit.Framework; using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using EngineLayer;
using MetaMorpheusCommandLine;

namespace Test
{
    [TestFixture]
    public static class TestCmd
    {
        [OneTimeSetUp]
        public static void Setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;
        }

        [Test]
        [Ignore("Ignored on AppVeyor")]
        public static void TestCommandLineMicroVignette()
        {
            //Stopwatch s = new Stopwatch();
            //s.Start();

            //string path = Path.Combine(TestContext.CurrentContext.TestDirectory, @"CommandLineMicroVignette");

            //// run the micro vignette via command-line
            //MetaMorpheusCommandLine.Program.Main(new string[] {
            //    "-v",
            //    "-o" + path } );

            //s.Stop();
            //Console.WriteLine("Command-line microvignette took: " + s.ToString());

            //Directory.Delete(path, true);
        }

        /// <summary>
        /// Regression test for nested .d folders (e.g., UnzippedFile.d/UnzippedFile.d).
        /// When users unzip timsTOF data, they sometimes end up with a nested structure
        /// where the outer .d folder is just a container and the inner .d folder is the actual data.
        /// </summary>
        [Test]
        public static void TestNestedDotDFolderDiscovery()
        {
            string tempFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestNestedDotD");

            try
            {
                // Clean up if exists from previous failed run
                if (Directory.Exists(tempFolder))
                {
                    Directory.Delete(tempFolder, true);
                }

                // Create nested .d structure: outer.d/inner.d with valid timsTOF files in inner
                string outerDotD = Path.Combine(tempFolder, "UnzippedFile.d");
                string innerDotD = Path.Combine(outerDotD, "UnzippedFile.d");
                Directory.CreateDirectory(innerDotD);

                // Create valid timsTOF files in the inner .d folder
                File.WriteAllText(Path.Combine(innerDotD, "analysis.tdf"), "fake tdf content");
                File.WriteAllText(Path.Combine(innerDotD, "analysis.tdf_bin"), "fake tdf_bin content");

                // Test the recursive discovery - should find the inner valid .d folder
                List<string> spectraFiles = new List<string>();

                // Use reflection or create a test-accessible method to call FindSpectraFilesRecursive
                // For now, we'll test the public behavior through CommandLineSettings
                var settings = new CommandLineSettings
                {
                    _spectra = new[] { tempFolder },
                    _tasks = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Task1-SearchTaskconfig.toml") },
                    _databases = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\gapdh.fasta") }
                };

                // This should not throw and should find the nested .d folder
                settings.ValidateCommandLineSettings();

                // Verify the inner .d folder was found
                Assert.That(settings.Spectra.Any(s => s.EndsWith("UnzippedFile.d")), Is.True,
                    "Should find the nested valid .d folder");

                // The outer invalid .d folder should not be in the list
                var dotDFolders = settings.Spectra.Where(s => s.EndsWith(".d")).ToList();
                Assert.That(dotDFolders.Count, Is.EqualTo(1), 
                    "Should find exactly one .d folder (the valid nested one)");
                Assert.That(dotDFolders[0], Does.EndWith(Path.Combine("UnzippedFile.d", "UnzippedFile.d")),
                    "The found .d folder should be the inner nested one");
            }
            finally
            {
                // Cleanup
                if (Directory.Exists(tempFolder))
                {
                    Directory.Delete(tempFolder, true);
                }
            }
        }

        /// <summary>
        /// Test that files inside a valid .d folder are not picked up separately.
        /// </summary>
        [Test]
        public static void TestFilesInsideValidDotDFolderAreNotAddedSeparately()
        {
            string tempFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestDotDFilesNotAddedCmd");

            try
            {
                // Clean up if exists from previous failed run
                if (Directory.Exists(tempFolder))
                {
                    Directory.Delete(tempFolder, true);
                }

                // Create a valid .d folder with additional files inside that should NOT be picked up
                string validDotD = Path.Combine(tempFolder, "valid_data.d");
                Directory.CreateDirectory(validDotD);

                // Create valid timsTOF files
                File.WriteAllText(Path.Combine(validDotD, "analysis.tdf"), "fake tdf content");
                File.WriteAllText(Path.Combine(validDotD, "analysis.tdf_bin"), "fake tdf_bin content");

                // Create other files inside that should NOT be added separately
                File.WriteAllText(Path.Combine(validDotD, "some_file.mzML"), "<mzML></mzML>");

                var settings = new CommandLineSettings
                {
                    _spectra = new[] { tempFolder },
                    _tasks = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Task1-SearchTaskconfig.toml") },
                    _databases = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\gapdh.fasta") }
                };

                settings.ValidateCommandLineSettings();

                // Should only have the .d folder, not the .mzML file inside it
                Assert.That(settings.Spectra.Count, Is.EqualTo(1), "Should find exactly one spectra entry");
                Assert.That(settings.Spectra[0], Does.EndWith("valid_data.d"), 
                    "Should be the .d folder, not files inside it");
                Assert.That(settings.Spectra.Any(s => s.EndsWith(".mzML")), Is.False,
                    "Should NOT find the mzML file inside the .d folder");
            }
            finally
            {
                // Cleanup
                if (Directory.Exists(tempFolder))
                {
                    Directory.Delete(tempFolder, true);
                }
            }
        }
    }
}
