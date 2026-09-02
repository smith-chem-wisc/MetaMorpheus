using EngineLayer;
using MetaMorpheusCommandLine;
using NUnit.Framework; using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Nett;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// The command line's design-file gate. These branches decide whether a search starts at all --
    /// and every one of them ends in either exit code 5 or a deleted file, so getting one wrong is
    /// either a run that refuses to start or a design file silently thrown away.
    /// </summary>
    [TestFixture]
    internal static class ResolveExperimentalDesignTests
    {
        private const int Continue = 0;
        private const int Stop = 5;

        private static string NewFolder(string name)
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, name);
            if (Directory.Exists(folder)) Directory.Delete(folder, true);
            Directory.CreateDirectory(folder);
            return folder;
        }

        /// <summary>A spectra file and a valid classic design naming it.</summary>
        private static string WriteClassicDesign(string folder, string spectraFile, bool valid)
        {
            string path = Path.Combine(folder, GlobalVariables.ExperimentalDesignFileName);
            File.WriteAllLines(path, valid
                ? new[] { "FileName\tCondition\tBiorep\tFraction\tTechrep",
                          $"{Path.GetFileName(spectraFile)}\tCondA\t1\t1\t1" }
                : new[] { "FileName\tCondition\tBiorep\tFraction\tTechrep",
                          $"{Path.GetFileName(spectraFile)}\tCondA\tnot-a-number\t1\t1" });
            return path;
        }

        private static string WriteTmtDesign(string folder, string spectraFile, bool valid)
        {
            string path = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);
            File.WriteAllLines(path, valid
                ? new[] { TmtExperimentalDesign.Header,
                          $"{spectraFile}\tPlexA\tS1\t126\tC1\t1\t1\t1" }
                : new[] { TmtExperimentalDesign.Header,
                          $"{spectraFile}\tPlexA\tS1\t126\tC1\tnot-a-number\t1\t1" });
            return path;
        }

        /// <summary>
        /// Two design files carrying replicate structure can drift apart, and nothing in the CLI could
        /// say which one the user meant, so the run stops rather than picking.
        /// </summary>
        [Test]
        public static void BothDesignFilesPresentIsRefused()
        {
            string folder = NewFolder("CmdDesignBoth");
            string spectra = Path.Combine(folder, "file.1.mzML");
            WriteClassicDesign(folder, spectra, valid: true);
            WriteTmtDesign(folder, spectra, valid: true);

            var written = new List<string>();
            int code = Program.ResolveExperimentalDesign(
                folder, new List<string> { spectra }, normalizationRequested: false,
                reportToConsole: true, write: written.Add, readLine: () => "n");

            Assert.That(code == Stop);
            Assert.That(written.Single().Contains("Only one design file type may be present"));

            Directory.Delete(folder, true);
        }

        /// <summary>Refused whether or not normalization was asked for -- the ambiguity is the problem.</summary>
        [Test]
        public static void BothDesignFilesPresentIsRefusedEvenWithoutNormalization()
        {
            string folder = NewFolder("CmdDesignBothQuiet");
            string spectra = Path.Combine(folder, "file.1.mzML");
            WriteClassicDesign(folder, spectra, valid: true);
            WriteTmtDesign(folder, spectra, valid: true);

            var written = new List<string>();
            int code = Program.ResolveExperimentalDesign(
                folder, new List<string> { spectra }, normalizationRequested: true,
                reportToConsole: false, write: written.Add, readLine: () => "n");

            Assert.That(code == Stop);
            Assert.That(written.Count == 0, "Nothing is printed at verbosity none");

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// A design is optional until something needs it. Normalization is that something: without one
        /// there is nothing to normalize across, so the run stops instead of quietly not normalizing.
        /// </summary>
        [Test]
        public static void NoDesignStopsOnlyWhenNormalizing()
        {
            string folder = NewFolder("CmdDesignNone");
            string spectra = Path.Combine(folder, "file.1.mzML");

            var written = new List<string>();
            Assert.That(Program.ResolveExperimentalDesign(
                folder, new List<string> { spectra }, normalizationRequested: false,
                reportToConsole: true, write: written.Add) == Continue,
                "No design and no normalization is a normal run");
            Assert.That(written.Count == 0);

            int code = Program.ResolveExperimentalDesign(
                folder, new List<string> { spectra }, normalizationRequested: true,
                reportToConsole: true, write: written.Add);

            Assert.That(code == Stop);
            Assert.That(written.Single().Contains("Normalization requires a design"));

            Directory.Delete(folder, true);
        }

        [Test]
        public static void NoDesignAndNormalizingIsSilentAtVerbosityNone()
        {
            string folder = NewFolder("CmdDesignNoneQuiet");
            var written = new List<string>();

            int code = Program.ResolveExperimentalDesign(
                folder, new List<string> { Path.Combine(folder, "file.1.mzML") },
                normalizationRequested: true, reportToConsole: false, write: written.Add);

            Assert.That(code == Stop);
            Assert.That(written.Count == 0);

            Directory.Delete(folder, true);
        }

        [TestCase(true, TestName = "ReadableClassicDesignIsAccepted")]
        [TestCase(false, TestName = "ReadableTmtDesignIsAccepted")]
        public static void AReadableDesignIsAcceptedAndNamed(bool classic)
        {
            string folder = NewFolder("CmdDesignGood" + (classic ? "Classic" : "Tmt"));
            string spectra = Path.Combine(folder, "file.1.mzML");
            if (classic) WriteClassicDesign(folder, spectra, valid: true);
            else WriteTmtDesign(folder, spectra, valid: true);

            var written = new List<string>();
            int code = Program.ResolveExperimentalDesign(
                folder, new List<string> { spectra }, normalizationRequested: true,
                reportToConsole: true, write: written.Add);

            Assert.That(code == Continue);
            Assert.That(written.Single() == "Read "
                + (classic ? GlobalVariables.ExperimentalDesignFileName : GlobalVariables.TmtExperimentalDesignFileName)
                + " successfully");

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// Normalizing against a design that could not be read would silently produce different
        /// numbers, so this is the one bad-design case that cannot be recovered from.
        /// </summary>
        [TestCase(true, TestName = "UnreadableClassicDesignStopsANormalizingRun")]
        [TestCase(false, TestName = "UnreadableTmtDesignStopsANormalizingRun")]
        public static void AnUnreadableDesignStopsANormalizingRunAndIsKept(bool classic)
        {
            string folder = NewFolder("CmdDesignBad" + (classic ? "Classic" : "Tmt"));
            string spectra = Path.Combine(folder, "file.1.mzML");
            string designPath = classic
                ? WriteClassicDesign(folder, spectra, valid: false)
                : WriteTmtDesign(folder, spectra, valid: false);

            var written = new List<string>();
            int code = Program.ResolveExperimentalDesign(
                folder, new List<string> { spectra }, normalizationRequested: true,
                reportToConsole: true, write: written.Add);

            Assert.That(code == Stop);
            Assert.That(written.Any(), "The errors are reported");
            Assert.That(File.Exists(designPath), "A design is never deleted out from under a normalizing run");

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// Not normalizing, so an unreadable design is recoverable -- by not reading it. Agreeing
        /// continues the run and LEAVES THE FILE ALONE. The prompt only ever asked about continuing,
        /// so deleting was an answer to a question nobody was asked; see #2256, and #2780 for the
        /// same change on the GUI side.
        /// </summary>
        [TestCase("y")]
        [TestCase("YES")]
        public static void AgreeingContinuesWithoutTheDesignAndKeepsTheFile(string answer)
        {
            string folder = NewFolder("CmdDesignContinue" + answer);
            string spectra = Path.Combine(folder, "file.1.mzML");
            string designPath = WriteClassicDesign(folder, spectra, valid: false);

            var written = new List<string>();
            int code = Program.ResolveExperimentalDesign(
                folder, new List<string> { spectra }, normalizationRequested: false,
                reportToConsole: true, write: written.Add, readLine: () => answer);

            Assert.That(code == Continue);
            Assert.That(File.Exists(designPath), "Continuing must never destroy the user's design");
            Assert.That(written.Any(w => w.Contains("Continue without an experimental design?")));
            Assert.That(written.Any(w => w.StartsWith("First error: ")));

            Directory.Delete(folder, true);
        }

        [TestCase("n")]
        [TestCase("")]
        [TestCase(null)]
        public static void DecliningStopsTheRunAndKeepsTheFile(string answer)
        {
            string folder = NewFolder("CmdDesignKeep" + (answer ?? "null").Length);
            string spectra = Path.Combine(folder, "file.1.mzML");
            string designPath = WriteClassicDesign(folder, spectra, valid: false);

            int code = Program.ResolveExperimentalDesign(
                folder, new List<string> { spectra }, normalizationRequested: false,
                reportToConsole: true, write: _ => { }, readLine: () => answer);

            Assert.That(code == Stop);
            Assert.That(File.Exists(designPath), "Declining keeps the design");

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// At verbosity "none" there is nobody to ask, so the recoverable case recovers itself -- by
        /// continuing without the design, NOT by deleting it. This previously destroyed the file with
        /// no prompt and no record, which meant `-v none` changed what the run did to the user's
        /// input rather than only what it printed. A verbosity flag does not carry that authority.
        /// </summary>
        [Test]
        public static void AtVerbosityNoneAnUnreadableDesignIsSkippedNotDeleted()
        {
            string folder = NewFolder("CmdDesignQuietSkip");
            string spectra = Path.Combine(folder, "file.1.mzML");
            string designPath = WriteTmtDesign(folder, spectra, valid: false);
            bool asked = false;

            int code = Program.ResolveExperimentalDesign(
                folder, new List<string> { spectra }, normalizationRequested: false,
                reportToConsole: false, write: _ => { }, readLine: () => { asked = true; return "n"; });

            Assert.That(code == Continue);
            Assert.That(!asked, "There is nobody to prompt at verbosity none");
            Assert.That(File.Exists(designPath), "-v none must not delete the user's design file");

            Directory.Delete(folder, true);
        }
        /// <summary>
        /// Every other test in this fixture calls ResolveExperimentalDesign directly, which leaves the
        /// wiring unproven: Run could drop the returned exit code, never call the gate, or resolve the
        /// design directory from somewhere other than the spectra. This drives
        /// MetaMorpheusCommandLine.Program.Main so the gate is reached the way a user reaches it.
        ///
        /// The design file here is malformed rather than absent, and that is what makes the assertion
        /// discriminating. An absent design stops the run from any directory -- so a gate pointed at
        /// the working directory would still return 5, and a test asserting only the exit code would
        /// pass while the directory logic was wrong. A malformed design beside the spectra separates
        /// them: read from the right directory the console carries that file's own parse error, and
        /// read from anywhere else it carries "no design present" instead.
        ///
        /// Exit code 5 and not 2 is what separates the design gate from settings validation, the other
        /// way this argument list can stop early. The refusal is also the only path safely drivable
        /// through Main, since every other outcome returns 0 and starts an actual search -- which is
        /// why the .mzML can be empty. That the run stops before anything opens it is the point.
        /// </summary>
        [Test]
        [NonParallelizable]
        public static void MainStopsANormalizingSearchWhoseDesignFileCannotBeRead()
        {
            string folder = NewFolder("CmdMainDesignGate");

            string spectra = Path.Combine(folder, "file.1.mzML");
            File.WriteAllText(spectra, string.Empty);

            string database = Path.Combine(folder, "db.fasta");
            File.WriteAllText(database, string.Empty);

            // Beside the spectra, and unreadable. Only a gate looking here reports this file's error.
            WriteTmtDesign(folder, spectra, valid: false);

            // Normalization is the one thing that makes a design file mandatory.
            SearchTask searchTask = new SearchTask();
            searchTask.SearchParameters.Normalize = true;
            string taskPath = Path.Combine(folder, "SearchTask.toml");
            Toml.WriteFile(searchTask, taskPath, MetaMorpheusTask.tomlConfig);

            string output = Path.Combine(folder, "output");

            TextWriter originalOut = Console.Out;
            StringWriter captured = new StringWriter();
            int exitCode;

            try
            {
                Console.SetOut(captured);
                exitCode = Program.Main(new[] { "-s", spectra, "-d", database, "-t", taskPath, "-o", output });
            }
            finally
            {
                Console.SetOut(originalOut);
            }

            Assert.That(exitCode == Stop,
                "A normalizing search whose design cannot be read must stop, and stop with the gate's code");
            Assert.That(captured.ToString().Contains("Biological Replicate"),
                "The design beside the spectra was not the one read -- the gate resolved a different directory");
            Assert.That(!Directory.Exists(output),
                "The search ran anyway -- Run did not honour the gate's exit code");

            Directory.Delete(folder, true);
        }
    }
}
