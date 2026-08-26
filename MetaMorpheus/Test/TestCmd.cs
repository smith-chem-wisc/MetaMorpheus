using NUnit.Framework; using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CommandLine;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using EngineLayer;
using MetaMorpheusCommandLine;
using Nett;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class TestCmd
    {
        private static GlobalSettings GlobalSettingsBeforeTest;
        private static string ScratchDataDirectory;

        [OneTimeSetUp]
        public static void Setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;
        }

        /// <summary>
        /// The licence tests below write a settings.toml and set GlobalVariables.GlobalSettings. Each gets
        /// a private scratch directory to stand in for the data directory, and the process-wide settings
        /// object is put back afterwards, so nothing they do is visible to the rest of the suite.
        /// </summary>
        [SetUp]
        public static void PerTestSetup()
        {
            GlobalSettingsBeforeTest = GlobalVariables.GlobalSettings;
            ScratchDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "ThermoLicenceScratch_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(ScratchDataDirectory);
        }

        [TearDown]
        public static void PerTestTearDown()
        {
            GlobalVariables.GlobalSettings = GlobalSettingsBeforeTest;

            if (Directory.Exists(ScratchDataDirectory))
            {
                Directory.Delete(ScratchDataDirectory, true);
            }
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

        /// <summary>
        /// Test that when a user provides a path to a .tdf file directly,
        /// the CMD corrects it to use the parent .d directory instead.
        /// </summary>
        [Test]
        public static void TestTdfFilePathCorrectedToParentDotDDirectory()
        {
            string tempFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestTdfPathCorrection");

            try
            {
                // Clean up if exists from previous failed run
                if (Directory.Exists(tempFolder))
                {
                    Directory.Delete(tempFolder, true);
                }

                // Create a valid .d folder with timsTOF files
                string validDotD = Path.Combine(tempFolder, "data.d");
                Directory.CreateDirectory(validDotD);

                string tdfFilePath = Path.Combine(validDotD, "analysis.tdf");
                string tdfBinFilePath = Path.Combine(validDotD, "analysis.tdf_bin");
                File.WriteAllText(tdfFilePath, "fake tdf content");
                File.WriteAllText(tdfBinFilePath, "fake tdf_bin content");

                // User provides the .tdf file path directly instead of the .d folder
                var settings = new CommandLineSettings
                {
                    _spectra = new[] { tdfFilePath },
                    _tasks = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Task1-SearchTaskconfig.toml") },
                    _databases = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\gapdh.fasta") }
                };

                settings.ValidateCommandLineSettings();

                // Should have corrected the path to the parent .d directory
                Assert.That(settings.Spectra.Count, Is.EqualTo(1), "Should have exactly one spectra entry");
                Assert.That(settings.Spectra[0], Does.EndWith("data.d"),
                    "Should have corrected the .tdf path to the parent .d directory");
                Assert.That(settings.Spectra[0], Does.Not.Contain("analysis.tdf"),
                    "Should not contain the .tdf filename");
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
        /// Test that when a user provides a path to a .tdf_bin file directly,
        /// the CMD corrects it to use the parent .d directory instead.
        /// </summary>
        [Test]
        public static void TestTdfBinFilePathCorrectedToParentDotDDirectory()
        {
            string tempFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestTdfBinPathCorrection");

            try
            {
                // Clean up if exists from previous failed run
                if (Directory.Exists(tempFolder))
                {
                    Directory.Delete(tempFolder, true);
                }

                // Create a valid .d folder with timsTOF files
                string validDotD = Path.Combine(tempFolder, "data.d");
                Directory.CreateDirectory(validDotD);

                string tdfFilePath = Path.Combine(validDotD, "analysis.tdf");
                string tdfBinFilePath = Path.Combine(validDotD, "analysis.tdf_bin");
                File.WriteAllText(tdfFilePath, "fake tdf content");
                File.WriteAllText(tdfBinFilePath, "fake tdf_bin content");

                // User provides the .tdf_bin file path directly instead of the .d folder
                var settings = new CommandLineSettings
                {
                    _spectra = new[] { tdfBinFilePath },
                    _tasks = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Task1-SearchTaskconfig.toml") },
                    _databases = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\gapdh.fasta") }
                };

                settings.ValidateCommandLineSettings();

                // Should have corrected the path to the parent .d directory
                Assert.That(settings.Spectra.Count, Is.EqualTo(1), "Should have exactly one spectra entry");
                Assert.That(settings.Spectra[0], Does.EndWith("data.d"),
                    "Should have corrected the .tdf_bin path to the parent .d directory");
                Assert.That(settings.Spectra[0], Does.Not.Contain("analysis.tdf_bin"),
                    "Should not contain the .tdf_bin filename");
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
        /// Regression test that a regular (non-tims) Bruker TOF .d folder, identified by an
        /// analysis.baf file, is discovered as a spectra file when its parent directory is provided.
        /// </summary>
        [Test]
        public static void TestBrukerBafDotDFolderDiscovery()
        {
            string tempFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestBrukerBafDiscovery");

            try
            {
                // Clean up if exists from previous failed run
                if (Directory.Exists(tempFolder))
                {
                    Directory.Delete(tempFolder, true);
                }

                // Create a valid regular-TOF .d folder (identified by analysis.baf)
                string validDotD = Path.Combine(tempFolder, "baf_data.d");
                Directory.CreateDirectory(validDotD);
                File.WriteAllText(Path.Combine(validDotD, "analysis.baf"), "fake baf content");

                // Also drop a file inside the .d that must NOT be picked up separately
                File.WriteAllText(Path.Combine(validDotD, "some_file.mzML"), "<mzML></mzML>");

                var settings = new CommandLineSettings
                {
                    _spectra = new[] { tempFolder },
                    _tasks = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Task1-SearchTaskconfig.toml") },
                    _databases = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\gapdh.fasta") }
                };

                settings.ValidateCommandLineSettings();

                Assert.That(settings.Spectra.Count, Is.EqualTo(1), "Should find exactly one spectra entry");
                Assert.That(settings.Spectra[0], Does.EndWith("baf_data.d"),
                    "Should discover the regular-TOF .d folder");
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

        /// <summary>
        /// Test that when a user provides a path to an analysis.baf file directly,
        /// the CMD corrects it to use the parent .d directory instead.
        /// </summary>
        [Test]
        public static void TestBafFilePathCorrectedToParentDotDDirectory()
        {
            string tempFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestBafPathCorrection");

            try
            {
                // Clean up if exists from previous failed run
                if (Directory.Exists(tempFolder))
                {
                    Directory.Delete(tempFolder, true);
                }

                // Create a valid regular-TOF .d folder (identified by analysis.baf)
                string validDotD = Path.Combine(tempFolder, "data.d");
                Directory.CreateDirectory(validDotD);

                string bafFilePath = Path.Combine(validDotD, "analysis.baf");
                File.WriteAllText(bafFilePath, "fake baf content");

                // User provides the .baf file path directly instead of the .d folder
                var settings = new CommandLineSettings
                {
                    _spectra = new[] { bafFilePath },
                    _tasks = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Task1-SearchTaskconfig.toml") },
                    _databases = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\gapdh.fasta") }
                };

                settings.ValidateCommandLineSettings();

                Assert.That(settings.Spectra.Count, Is.EqualTo(1), "Should have exactly one spectra entry");
                Assert.That(settings.Spectra[0], Does.EndWith("data.d"),
                    "Should have corrected the .baf path to the parent .d directory");
                Assert.That(settings.Spectra[0], Does.Not.Contain("analysis.baf"),
                    "Should not contain the .baf filename");
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
        /// Test that Bruker inner files encountered during recursive directory search are collapsed
        /// to their parent .d folder. A trailing directory separator on the .d path (which shells add
        /// during tab-completion) prevents the .d folder from being recognized up front, so the search
        /// recurses into it and reaches analysis.tdf/analysis.tdf_bin as individual files.
        /// Both inner files map to the same .d folder, so the result must also be de-duplicated.
        /// </summary>
        [Test]
        public static void TestInnerBrukerFilesFoundByRecursionCollapseToParentDotDDirectory()
        {
            string tempFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestRecursiveInnerFileCollapse");

            try
            {
                // Clean up if exists from previous failed run
                if (Directory.Exists(tempFolder))
                {
                    Directory.Delete(tempFolder, true);
                }

                // Create a valid timsTOF .d folder
                string validDotD = Path.Combine(tempFolder, "data.d");
                Directory.CreateDirectory(validDotD);
                File.WriteAllText(Path.Combine(validDotD, "analysis.tdf"), "fake tdf content");
                File.WriteAllText(Path.Combine(validDotD, "analysis.tdf_bin"), "fake tdf_bin content");

                // The .d folder is given with a trailing separator, e.g. "data.d\"
                var settings = new CommandLineSettings
                {
                    _spectra = new[] { validDotD + Path.DirectorySeparatorChar },
                    _tasks = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Task1-SearchTaskconfig.toml") },
                    _databases = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\gapdh.fasta") },
                    Verbosity = CommandLineSettings.VerbosityType.normal
                };

                settings.ValidateCommandLineSettings();

                Assert.That(settings.Spectra.Count, Is.EqualTo(1),
                    "Both inner files should collapse to a single .d folder entry");
                Assert.That(settings.Spectra[0], Does.EndWith("data.d"),
                    "Should have collapsed the inner Bruker files to the parent .d directory");
                Assert.That(settings.Spectra.Any(s => s.Contains("analysis.tdf")), Is.False,
                    "Should not contain the individual inner Bruker files");
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
        /// Test that a Bruker inner file inside an *invalid* .d folder is not promoted to that folder.
        /// An incomplete timsTOF folder (analysis.tdf with no analysis.tdf_bin) is not a readable .d
        /// folder, so the .tdf file itself is kept. Non-spectra files in the same folder are ignored.
        /// </summary>
        [Test]
        public static void TestInnerBrukerFileInInvalidDotDFolderIsNotCollapsed()
        {
            string tempFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestInvalidDotDInnerFile");

            try
            {
                // Clean up if exists from previous failed run
                if (Directory.Exists(tempFolder))
                {
                    Directory.Delete(tempFolder, true);
                }

                // Create an incomplete .d folder: analysis.tdf present, analysis.tdf_bin missing
                string invalidDotD = Path.Combine(tempFolder, "incomplete.d");
                Directory.CreateDirectory(invalidDotD);
                File.WriteAllText(Path.Combine(invalidDotD, "analysis.tdf"), "fake tdf content");

                // A non-spectra file that must be ignored by the recursive search
                File.WriteAllText(Path.Combine(invalidDotD, "notes.txt"), "not a spectra file");

                var settings = new CommandLineSettings
                {
                    _spectra = new[] { tempFolder },
                    _tasks = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Task1-SearchTaskconfig.toml") },
                    _databases = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\gapdh.fasta") }
                };

                settings.ValidateCommandLineSettings();

                Assert.That(settings.Spectra.Count, Is.EqualTo(1), "Should find exactly one spectra entry");
                Assert.That(settings.Spectra[0], Does.EndWith(Path.Combine("incomplete.d", "analysis.tdf")),
                    "Should keep the .tdf file, since its parent .d folder is not a readable Bruker folder");
                Assert.That(settings.Spectra.Any(s => s.EndsWith(".txt")), Is.False,
                    "Should not pick up non-spectra files");
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
        /// -g must emit a default toml for every MyTask type. Averaging was omitted, which left
        /// spectral averaging as the one task whose settings could not be discovered from the
        /// command line at all - the enum has six members and Program.cs runs all six.
        /// </summary>
        [Test]
        public static void TestGenerateDefaultTaskTomlsCoversEveryTaskType()
        {
            string tempFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestGenerateDefaultTomls");

            try
            {
                if (Directory.Exists(tempFolder))
                {
                    Directory.Delete(tempFolder, true);
                }

                CommandLineSettings.GenerateDefaultTaskTomls(tempFolder);

                string[] expected =
                {
                    "CalibrationTask.toml", "GptmdTask.toml", "SearchTask.toml",
                    "XLSearchTask.toml", "GlycoSearchTask.toml", "AveragingTask.toml"
                };

                foreach (string name in expected)
                {
                    Assert.That(File.Exists(Path.Combine(tempFolder, name)), Is.True,
                        $"-g did not write {name}");
                }

                // One file per MyTask member, no more: a task type added to the enum without a
                // generated default is exactly the gap this test exists to catch.
                Assert.That(Directory.GetFiles(tempFolder, "*.toml").Length,
                    Is.EqualTo(Enum.GetValues(typeof(MyTask)).Length),
                    "-g should write exactly one default toml per MyTask member");
            }
            finally
            {
                if (Directory.Exists(tempFolder))
                {
                    Directory.Delete(tempFolder, true);
                }
            }
        }

        /// <summary>
        /// The generated averaging toml must be runnable, not merely present. Unlike its five
        /// siblings, SpectralAveragingTask's parameterless constructor initialised nothing, so a
        /// default-constructed task would have serialised null sections - a file that exists,
        /// looks plausible, and fails when someone tries to use it.
        /// </summary>
        [Test]
        public static void TestGeneratedAveragingTomlRoundTripsWithPopulatedParameters()
        {
            string tempFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestAveragingTomlRoundTrip");

            try
            {
                if (Directory.Exists(tempFolder))
                {
                    Directory.Delete(tempFolder, true);
                }

                CommandLineSettings.GenerateDefaultTaskTomls(tempFolder);
                string averagingToml = Path.Combine(tempFolder, "AveragingTask.toml");

                // The TaskType string is what Program.cs switches on to dispatch the task.
                var raw = Toml.ReadFile(averagingToml, MetaMorpheusTask.tomlConfig);
                Assert.That(raw.Get<string>("TaskType"), Is.EqualTo(MyTask.Average.ToString()));

                var task = Toml.ReadFile<SpectralAveragingTask>(averagingToml, MetaMorpheusTask.tomlConfig);
                Assert.That(task.Parameters, Is.Not.Null, "generated toml carried no [Parameters]");
                Assert.That(task.CommonParameters, Is.Not.Null, "generated toml carried no [CommonParameters]");
                Assert.That(task.TaskType, Is.EqualTo(MyTask.Average));
            }
            finally
            {
                if (Directory.Exists(tempFolder))
                {
                    Directory.Delete(tempFolder, true);
                }
            }
        }

        /// <summary>
        /// The option name is the whole interface here - it is what goes in a Dockerfile and a job
        /// script - so a rename or a typo in the attribute has to fail loudly rather than quietly
        /// leaving the flag unbound and the run back at the interactive prompt.
        /// </summary>
        [Test]
        public static void TestAcceptThermoLicenceFlagBindsFromTheCommandLine()
        {
            bool bound = false;

            CommandLine.Parser.Default
                .ParseArguments<CommandLineSettings>(new[] { "--acceptThermoLicence" })
                .WithParsed(options => bound = options.AcceptThermoLicence);

            Assert.That(bound, Is.True, "--acceptThermoLicence was not bound by the command line parser");
        }

        /// <summary>
        /// --acceptThermoLicence on its own is a setup step, not a run, so it must not be held to the
        /// task/database/spectra requirements. Without the exemption the invocation that a container
        /// build or a first-time macOS user needs most is rejected before it reaches Program.Run.
        /// </summary>
        [Test]
        public static void TestAcceptThermoLicenceValidatesWithoutARunSpecified()
        {
            var settings = new CommandLineSettings { AcceptThermoLicence = true };

            Assert.DoesNotThrow(() => settings.ValidateCommandLineSettings());

            Assert.That(settings.Tasks.Count, Is.EqualTo(0));
            Assert.That(settings.Databases.Count, Is.EqualTo(0));
            Assert.That(settings.Spectra.Count, Is.EqualTo(0));
        }

        /// <summary>
        /// The exemption above is conditional on the flag. An empty command line is still an error.
        /// </summary>
        [Test]
        public static void TestEmptyCommandLineIsStillRejectedWithoutTheLicenceFlag()
        {
            var settings = new CommandLineSettings();

            Assert.Throws<MetaMorpheusException>(() => settings.ValidateCommandLineSettings());
        }

        /// <summary>
        /// Given alongside a run the flag is only a modifier, so the run itself is validated as usual.
        /// A half-specified search must not become valid just because the licence was accepted with it.
        /// </summary>
        [Test]
        public static void TestAcceptThermoLicenceDoesNotExemptARunFromValidation()
        {
            var settings = new CommandLineSettings
            {
                AcceptThermoLicence = true,
                _databases = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\gapdh.fasta") }
            };

            var exception = Assert.Throws<MetaMorpheusException>(() => settings.ValidateCommandLineSettings());

            Assert.That(exception.Message, Does.Contain("task"),
                "a run given with --acceptThermoLicence should still be rejected for having no task");
        }

        /// <summary>
        /// Without the flag this must do nothing at all - not print, not write, and above all not stop
        /// the run - since it sits on the path of every single command line invocation.
        /// </summary>
        [Test]
        public static void TestAgreeToThermoLicenceIsANoOpWithoutTheFlag()
        {
            var settings = new CommandLineSettings { AcceptThermoLicence = false };
            var output = new StringWriter();

            bool nothingToRun = Program.AgreeToThermoLicence(settings, ScratchDataDirectory, output);

            Assert.That(nothingToRun, Is.False, "a run without the flag must not be stopped");
            Assert.That(output.ToString(), Is.Empty);
            Assert.That(File.Exists(Path.Combine(ScratchDataDirectory, "settings.toml")), Is.False,
                "nothing should have been written without the flag");
        }

        /// <summary>
        /// The flag on an unagreed machine: the licence is printed, the agreement is recorded on disk and
        /// applied to the running process, and the caller is told there is no run to continue into.
        /// WriteExcelCompatibleTSVs starts true and is asserted true afterwards because it is the other
        /// half of settings.toml - recording an agreement must not quietly reset an unrelated setting.
        /// </summary>
        [Test]
        public static void TestAgreeToThermoLicencePrintsAndRecordsWhenNotYetAgreed()
        {
            GlobalVariables.GlobalSettings = new GlobalSettings
            {
                UserHasAgreedToThermoRawFileReaderLicence = false,
                WriteExcelCompatibleTSVs = true
            };

            var settings = new CommandLineSettings
            {
                AcceptThermoLicence = true,
                Verbosity = CommandLineSettings.VerbosityType.normal
            };
            settings.ValidateCommandLineSettings();
            var output = new StringWriter();

            bool nothingToRun = Program.AgreeToThermoLicence(settings, ScratchDataDirectory, output);

            Assert.That(nothingToRun, Is.True, "the flag on its own is a setup step, not a run");
            Assert.That(output.ToString(), Does.Contain("SOFTWARE LICENSE AGREEMENT"),
                "the licence itself must be printed, so what was agreed to is in the log");
            Assert.That(output.ToString(), Does.Contain("Recorded in"));

            var recorded = Toml.ReadFile<GlobalSettings>(Path.Combine(ScratchDataDirectory, "settings.toml"));
            Assert.That(recorded.UserHasAgreedToThermoRawFileReaderLicence, Is.True,
                "the agreement was not written to settings.toml");
            Assert.That(recorded.WriteExcelCompatibleTSVs, Is.True,
                "recording the agreement reset the other setting in settings.toml");
            Assert.That(GlobalVariables.GlobalSettings.UserHasAgreedToThermoRawFileReaderLicence, Is.True,
                "the agreement was written but not applied to the running process");
        }

        /// <summary>
        /// On a machine that has already agreed there is nothing to record and nothing to re-read, so the
        /// licence is not printed again and settings.toml is left alone. It still reports success.
        /// </summary>
        [Test]
        public static void TestAgreeToThermoLicenceDoesNotReprintOnceAlreadyAgreed()
        {
            GlobalVariables.GlobalSettings = new GlobalSettings
            {
                UserHasAgreedToThermoRawFileReaderLicence = true,
                WriteExcelCompatibleTSVs = false
            };

            var settings = new CommandLineSettings
            {
                AcceptThermoLicence = true,
                Verbosity = CommandLineSettings.VerbosityType.normal
            };
            settings.ValidateCommandLineSettings();
            var output = new StringWriter();

            bool nothingToRun = Program.AgreeToThermoLicence(settings, ScratchDataDirectory, output);

            Assert.That(nothingToRun, Is.True);
            Assert.That(output.ToString(), Does.Not.Contain("SOFTWARE LICENSE AGREEMENT"),
                "the licence should not be reprinted to someone who has already agreed");
            Assert.That(output.ToString(), Does.Contain("Agreed to the Thermo RawFileReader licence"));
            Assert.That(File.Exists(Path.Combine(ScratchDataDirectory, "settings.toml")), Is.False,
                "an already-agreed machine should not have its settings.toml rewritten");
        }

        /// <summary>
        /// -v none means no output, and that applies to the confirmation as much as to anything else.
        /// </summary>
        [Test]
        public static void TestAgreeToThermoLicenceSaysNothingAtVerbosityNone()
        {
            GlobalVariables.GlobalSettings = new GlobalSettings { UserHasAgreedToThermoRawFileReaderLicence = true };

            var settings = new CommandLineSettings
            {
                AcceptThermoLicence = true,
                Verbosity = CommandLineSettings.VerbosityType.none
            };
            settings.ValidateCommandLineSettings();
            var output = new StringWriter();

            Program.AgreeToThermoLicence(settings, ScratchDataDirectory, output);

            Assert.That(output.ToString(), Is.Empty, "-v none should suppress the confirmation");
        }

        /// <summary>
        /// -v minimal is the other half of that condition, and does report the agreement.
        /// </summary>
        [Test]
        public static void TestAgreeToThermoLicenceReportsAtVerbosityMinimal()
        {
            GlobalVariables.GlobalSettings = new GlobalSettings { UserHasAgreedToThermoRawFileReaderLicence = true };

            var settings = new CommandLineSettings
            {
                AcceptThermoLicence = true,
                Verbosity = CommandLineSettings.VerbosityType.minimal
            };
            settings.ValidateCommandLineSettings();
            var output = new StringWriter();

            Program.AgreeToThermoLicence(settings, ScratchDataDirectory, output);

            Assert.That(output.ToString(), Does.Contain("Agreed to the Thermo RawFileReader licence"));
        }

        /// <summary>
        /// Given alongside a real run the flag is only a modifier: the agreement is recorded and the
        /// caller is told to carry on into the search rather than stopping at the setup step.
        /// </summary>
        [Test]
        public static void TestAgreeToThermoLicenceContinuesIntoARunWhenOneWasGiven()
        {
            GlobalVariables.GlobalSettings = new GlobalSettings { UserHasAgreedToThermoRawFileReaderLicence = false };

            var settings = new CommandLineSettings
            {
                AcceptThermoLicence = true,
                _tasks = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Task1-SearchTaskconfig.toml") },
                _databases = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\gapdh.fasta") },
                _spectra = new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch\SmallCalibratible_Yeast.mzML") },
                OutputFolder = ScratchDataDirectory
            };
            settings.ValidateCommandLineSettings();
            var output = new StringWriter();

            bool nothingToRun = Program.AgreeToThermoLicence(settings, ScratchDataDirectory, output);

            Assert.That(nothingToRun, Is.False, "a run was given, so the flag must not stop it");
            Assert.That(Toml.ReadFile<GlobalSettings>(Path.Combine(ScratchDataDirectory, "settings.toml"))
                .UserHasAgreedToThermoRawFileReaderLicence, Is.True,
                "the agreement should still have been recorded on the way through");
        }

        /// <summary>
        /// The point of the flag is what it writes, so this drives the real entry point rather than a
        /// stand-in: an unagreed settings.toml goes in, MetaMorpheusCommandLine.Program.Main runs with
        /// nothing but the flag, and the agreement comes out recorded on disk and in memory.
        ///
        /// Given on its own the flag returns before Run() subscribes any engine or task handlers, so the
        /// only global state this touches is settings.toml and GlobalVariables.GlobalSettings, both of
        /// which are restored below. WriteExcelCompatibleTSVs is set beforehand and asserted afterwards
        /// because it is the other half of that file: recording an agreement must not quietly reset the
        /// one unrelated setting sharing it.
        /// </summary>
        [Test]
        [NonParallelizable]
        public static void TestAcceptThermoLicenceRecordsTheAgreementThroughMain()
        {
            string settingsPath = Path.Combine(GlobalVariables.DataDir, @"settings.toml");
            string settingsBackup = File.Exists(settingsPath) ? File.ReadAllText(settingsPath) : null;
            GlobalSettings globalSettingsBackup = GlobalVariables.GlobalSettings;

            try
            {
                // a machine that has never agreed, and that has something to lose in the same file
                Toml.WriteFile(new GlobalSettings
                {
                    UserHasAgreedToThermoRawFileReaderLicence = false,
                    WriteExcelCompatibleTSVs = true
                }, settingsPath);

                int exitCode = Program.Main(new[] { "--acceptThermoLicence" });

                Assert.That(exitCode, Is.EqualTo(0), "--acceptThermoLicence on its own should succeed");

                var recorded = Toml.ReadFile<GlobalSettings>(settingsPath);
                Assert.That(recorded.UserHasAgreedToThermoRawFileReaderLicence, Is.True,
                    "the agreement was not written to settings.toml");
                Assert.That(recorded.WriteExcelCompatibleTSVs, Is.True,
                    "recording the agreement reset the other setting in settings.toml");
                Assert.That(GlobalVariables.GlobalSettings.UserHasAgreedToThermoRawFileReaderLicence, Is.True,
                    "the agreement was written but not applied to the running process");

                // Running it again on an already-agreed machine takes the other branch: nothing to
                // record, nothing to re-print, and still a success rather than a complaint.
                int secondExitCode = Program.Main(new[] { "--acceptThermoLicence" });

                Assert.That(secondExitCode, Is.EqualTo(0), "--acceptThermoLicence should be idempotent");
                Assert.That(Toml.ReadFile<GlobalSettings>(settingsPath).UserHasAgreedToThermoRawFileReaderLicence,
                    Is.True);
            }
            finally
            {
                if (settingsBackup == null)
                {
                    File.Delete(settingsPath);
                }
                else
                {
                    File.WriteAllText(settingsPath, settingsBackup);
                }

                GlobalVariables.GlobalSettings = globalSettingsBackup;
            }
        }
        /// <summary>
        /// Startup warnings collected during GlobalVariables.SetUpGlobalVariables (e.g. a failure to seed
        /// MonosaccharidesCustom.tsv) have to reach the console on a command-line run -- the GUI shows them in the
        /// notifications pane, and the CLI has nowhere else to put them. Written once, then cleared.
        /// </summary>
        [Test]
        [TestCase(CommandLineSettings.VerbosityType.minimal)]
        [TestCase(CommandLineSettings.VerbosityType.normal)]
        [NonParallelizable] // mutates the process-wide GlobalVariables.StartupWarnings
        public static void TestFlushStartupWarningsWritesAndClears(CommandLineSettings.VerbosityType verbosity)
        {
            var original = GlobalVariables.StartupWarnings.ToList();
            var originalOut = Console.Out;
            var sw = new StringWriter();

            try
            {
                GlobalVariables.StartupWarnings.Clear();
                GlobalVariables.StartupWarnings.Add("could not create the custom monosaccharide file");
                GlobalVariables.StartupWarnings.Add("second warning");
                Console.SetOut(sw);

                Program.FlushStartupWarnings(verbosity);
            }
            finally
            {
                Console.SetOut(originalOut);
                GlobalVariables.StartupWarnings.Clear();
                GlobalVariables.StartupWarnings.AddRange(original);
            }

            string written = sw.ToString();
            Assert.That(written, Does.Contain("could not create the custom monosaccharide file"));
            Assert.That(written, Does.Contain("second warning"));
        }

        /// <summary>
        /// -v none means no output at all, so the warnings are dropped rather than written. They are still cleared,
        /// so a later caller cannot report them a second time.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.StartupWarnings
        public static void TestFlushStartupWarningsSilentAtVerbosityNone()
        {
            var original = GlobalVariables.StartupWarnings.ToList();
            var originalOut = Console.Out;
            var sw = new StringWriter();
            int remaining;

            try
            {
                GlobalVariables.StartupWarnings.Clear();
                GlobalVariables.StartupWarnings.Add("could not create the custom monosaccharide file");
                Console.SetOut(sw);

                Program.FlushStartupWarnings(CommandLineSettings.VerbosityType.none);
                remaining = GlobalVariables.StartupWarnings.Count;
            }
            finally
            {
                Console.SetOut(originalOut);
                GlobalVariables.StartupWarnings.Clear();
                GlobalVariables.StartupWarnings.AddRange(original);
            }

            Assert.That(sw.ToString(), Is.Empty, "nothing should be written at verbosity 'none'");
            Assert.That(remaining, Is.EqualTo(0), "warnings should be cleared even when they are not written");
        }

        /// <summary>
        /// The common case: nothing went wrong at startup, so nothing is written.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.StartupWarnings
        public static void TestFlushStartupWarningsWritesNothingWhenThereAreNoWarnings()
        {
            var original = GlobalVariables.StartupWarnings.ToList();
            var originalOut = Console.Out;
            var sw = new StringWriter();

            try
            {
                GlobalVariables.StartupWarnings.Clear();
                Console.SetOut(sw);

                Program.FlushStartupWarnings(CommandLineSettings.VerbosityType.normal);
            }
            finally
            {
                Console.SetOut(originalOut);
                GlobalVariables.StartupWarnings.Clear();
                GlobalVariables.StartupWarnings.AddRange(original);
            }

            Assert.That(sw.ToString(), Is.Empty);
        }
    }
}
