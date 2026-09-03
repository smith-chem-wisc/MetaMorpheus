using NUnit.Framework; using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System;
using System.Diagnostics;
using System.IO;
using EngineLayer;

namespace Test
{
    /// <summary>
    /// The command line asks two yes/no questions, and both used to dereference
    /// Console.ReadLine() unguarded. When there is no console to answer from - a CI job, a
    /// container, a cluster batch job, anything with stdin redirected or closed - ReadLine()
    /// returns null and the run died with a NullReferenceException instead of explaining itself.
    ///
    /// These run the real CLI as a separate process rather than calling Program.Main in-process,
    /// for two reasons. It is the honest version of the condition under test: the process really
    /// has no stdin, rather than a StringReader standing in for one, and the assertion is on the
    /// process exit code a script would actually see. And Program.Run subscribes static handlers
    /// to MetaMorpheusEngine's events and never unsubscribes them, so calling it in-process leaves
    /// those handlers attached for the rest of the run - which made eleven unrelated parsimony
    /// tests fail later in the suite while passing on their own.
    /// </summary>
    [TestFixture]
    public class CmdConsolePromptTests
    {
        private string _settingsPath;
        private string _savedSettings;

        /// <summary>
        /// For a non-installed build GlobalVariables.DataDir is the build output directory, so
        /// settings.toml sits beside the binaries and the test owns it. The licence prompt only
        /// appears while it records no agreement, so each test starts from "not yet agreed".
        /// </summary>
        [SetUp]
        public void SetUp()
        {
            _settingsPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "settings.toml");
            _savedSettings = File.Exists(_settingsPath) ? File.ReadAllText(_settingsPath) : null;

            File.WriteAllText(_settingsPath,
                "UserHasAgreedToThermoRawFileReaderLicence = false" + Environment.NewLine
                + "WriteExcelCompatibleTSVs = true" + Environment.NewLine);
        }

        [TearDown]
        public void TearDown()
        {
            if (_savedSettings != null)
            {
                File.WriteAllText(_settingsPath, _savedSettings);
            }
            else if (File.Exists(_settingsPath))
            {
                File.Delete(_settingsPath);
            }
        }

        /// <summary>
        /// Runs the built CLI. A null answer closes stdin outright, which is what a redirected or
        /// closed stdin does to Console.ReadLine() in a real batch run.
        /// </summary>
        private static (int ExitCode, string Output) RunCli(string answer, params string[] args)
        {
            string testDirectory = TestContext.CurrentContext.TestDirectory;

            var startInfo = new ProcessStartInfo
            {
                FileName = "dotnet",
                WorkingDirectory = testDirectory,
                RedirectStandardInput = true,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                UseShellExecute = false,
                CreateNoWindow = true
            };

            startInfo.ArgumentList.Add(Path.Combine(testDirectory, "CMD.dll"));
            foreach (string arg in args)
            {
                startInfo.ArgumentList.Add(arg);
            }

            using Process process = Process.Start(startInfo);

            if (answer != null)
            {
                process.StandardInput.Write(answer);
            }
            process.StandardInput.Close();

            string output = process.StandardOutput.ReadToEnd() + process.StandardError.ReadToEnd();
            Assert.IsTrue(process.WaitForExit(300_000), "the CLI did not exit; it is most likely blocked on a prompt");

            return (process.ExitCode, output);
        }

        private static string MakeEmptyFile(string folder, string fileName)
        {
            string path = Path.Combine(folder, fileName);
            File.WriteAllBytes(path, Array.Empty<byte>());
            return path;
        }

        private static string TestDataPath(string fileName) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", fileName);

        private static string SearchTaskTomlPath() =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "SlicedSearchTaskConfig.toml");

        private static string FreshFolder(string name)
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, name);
            if (Directory.Exists(folder))
            {
                Directory.Delete(folder, true);
            }
            Directory.CreateDirectory(folder);
            return folder;
        }

        /// <summary>
        /// A .raw file on a machine that has not agreed to the Thermo licence, with nothing to read
        /// the answer from. Exit code 6 - NOT 5, which already means the experimental design failed
        /// at three separate sites, and NOT 3, which means the user actively declined.
        /// </summary>
        [Test]
        public void ThermoLicencePromptWithNoConsoleReturnsSix()
        {
            string folder = FreshFolder("CmdThermoLicenceNoConsole");

            try
            {
                // the file only has to exist and end in .raw; the licence gate is reached before
                // anything tries to read it
                string rawFile = MakeEmptyFile(folder, "doesNotMatter.raw");

                var result = RunCli(null,
                    "-t", SearchTaskTomlPath(),
                    "-s", rawFile,
                    "-d", TestDataPath("smalldb.fasta"),
                    "-o", Path.Combine(folder, "output"));

                Assert.AreEqual(6, result.ExitCode);
                Assert.IsTrue(result.Output.Contains("not running interactively"),
                    "the run should say why it could not ask, rather than failing silently");
                Assert.IsTrue(result.Output.Contains("UserHasAgreedToThermoRawFileReaderLicence"),
                    "the run should name the setting that lets a batch job agree without a console");
                Assert.IsFalse(result.Output.Contains("NullReferenceException"));

                // agreeing is recorded in settings.toml; refusing to answer must record nothing
                Assert.IsTrue(File.ReadAllText(_settingsPath).Contains("false"),
                    "an unanswerable prompt must not be treated as agreement");
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// Exit code 5 already means the experimental design failed. This pins that a licence prompt
        /// nobody could answer stays distinguishable from an active refusal - the whole reason for
        /// not reusing an existing code.
        /// </summary>
        [Test]
        public void DeclinedThermoLicenceStillReturnsThreeNotSix()
        {
            string folder = FreshFolder("CmdThermoLicenceDeclined");

            try
            {
                string rawFile = MakeEmptyFile(folder, "doesNotMatter.raw");

                var result = RunCli("n" + Environment.NewLine,
                    "-t", SearchTaskTomlPath(),
                    "-s", rawFile,
                    "-d", TestDataPath("smalldb.fasta"),
                    "-o", Path.Combine(folder, "output"));

                Assert.AreEqual(3, result.ExitCode, "a refusal is 3; only an unanswerable prompt is 6");
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// A piped answer arriving with trailing whitespace or a stray carriage return is the normal
        /// shape of `echo n | metamorpheus`, and used to be read as neither y nor n.
        /// </summary>
        [Test]
        public void PipedAnswerWithTrailingWhitespaceIsStillRead()
        {
            string folder = FreshFolder("CmdThermoLicenceWhitespace");

            try
            {
                string rawFile = MakeEmptyFile(folder, "doesNotMatter.raw");

                var result = RunCli("  n\t \r\n",
                    "-t", SearchTaskTomlPath(),
                    "-s", rawFile,
                    "-d", TestDataPath("smalldb.fasta"),
                    "-o", Path.Combine(folder, "output"));

                Assert.AreEqual(3, result.ExitCode, "a padded 'n' is a refusal, not an unreadable answer");
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// The second prompt: an experimental design that will not parse, and no console to ask
        /// whether it may be deleted. This is the site that still threw after the licence prompt was
        /// fixed. Refusing to answer must leave the user's file on disk.
        /// </summary>
        [Test]
        public void UnparsableExperimentalDesignWithNoConsoleReturnsFiveAndKeepsTheFile()
        {
            string folder = FreshFolder("CmdExperimentalDesignNoConsole");

            try
            {
                // .mzML, so the Thermo licence gate is not in the way of the prompt under test
                string spectraFile = MakeEmptyFile(folder, "spectra.mzML");

                // a data row with three cells where five are required
                string designPath = Path.Combine(folder, GlobalVariables.ExperimentalDesignFileName);
                File.WriteAllText(designPath,
                    "FileName\tCondition\tBiorep\tFraction\tTechrep" + Environment.NewLine
                    + "spectra.mzML\tcondition\t1" + Environment.NewLine);

                var result = RunCli(null,
                    "-t", SearchTaskTomlPath(),
                    "-s", spectraFile,
                    "-d", TestDataPath("smalldb.fasta"),
                    "-o", Path.Combine(folder, "output"));

                Assert.AreEqual(5, result.ExitCode);
                Assert.IsFalse(result.Output.Contains("NullReferenceException"),
                    "the unanswerable prompt must be reported, not thrown");
                Assert.IsTrue(File.Exists(designPath),
                    "an unanswerable prompt must not delete the user's experimental design");
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }
    }
}
