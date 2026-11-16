using EngineLayer;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test
{
    /// <summary>
    /// The design-file rule the command line and the GUI's Run button share. The command line's own output is pinned by
    /// ResolveExperimentalDesignTests, which pass unchanged; these pin the decision itself, above all the cases the GUI
    /// used to decide differently: it read only ExperimentalDesign.tsv, so it never refused both files, refused a
    /// normalizing TMT run whose design was TmtDesign.txt, and never noticed an unreadable TmtDesign.txt.
    /// </summary>
    [TestFixture]
    internal static class ExperimentalDesignCheckTests
    {
        private static string NewFolder(string name)
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "ExperimentalDesignCheck", name);
            if (Directory.Exists(folder)) Directory.Delete(folder, true);
            Directory.CreateDirectory(folder);
            return folder;
        }

        private static void WriteClassicDesign(string folder, string spectraFile, bool valid, params int[] bioreps)
        {
            if (bioreps.Length == 0) bioreps = new[] { 1 };
            File.WriteAllLines(Path.Combine(folder, GlobalVariables.ExperimentalDesignFileName),
                new[] { "FileName\tCondition\tBiorep\tFraction\tTechrep" }.Concat(bioreps.Select((b, i) =>
                    $"{Path.GetFileNameWithoutExtension(spectraFile)}{(i == 0 ? "" : "_" + i)}.mzML\tCondA\t{(valid ? b.ToString() : "not-a-number")}\t1\t1")));
        }

        private static void WriteTmtDesign(string folder, string spectraFile, bool valid)
        {
            File.WriteAllLines(Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName), new[]
            {
                TmtExperimentalDesign.Header,
                $"{spectraFile}\tPlexA\tS1\t126\tC1\t{(valid ? "1" : "not-a-number")}\t1\t1",
            });
        }

        [TestCase(false)]
        [TestCase(true)]
        public static void BothDesignFilesAreRefusedWhetherOrNotNormalizing(bool normalize)
        {
            string folder = NewFolder("Both" + normalize);
            string spectra = Path.Combine(folder, "file.1.mzML");
            WriteClassicDesign(folder, spectra, valid: true);
            WriteTmtDesign(folder, spectra, valid: true);

            var check = ExperimentalDesignCheck.Check(folder, new List<string> { spectra }, normalize);

            Assert.That(check.Verdict, Is.EqualTo(ExperimentalDesignVerdict.Refuse));
            Assert.That(check.Problem, Is.EqualTo(ExperimentalDesignProblem.BothDesignFiles));
            Assert.That(check.Messages.Single(), Does.Contain("Only one design file type may be present"));
        }

        [Test]
        public static void NoDesignIsRefusedOnlyWhenNormalizing()
        {
            string folder = NewFolder("None");
            var spectra = new List<string> { Path.Combine(folder, "file.1.mzML") };

            var quiet = ExperimentalDesignCheck.Check(folder, spectra, normalizationRequested: false);
            Assert.That((quiet.Verdict, quiet.Problem, quiet.DesignFileName), Is.EqualTo((ExperimentalDesignVerdict.Proceed, ExperimentalDesignProblem.None, (string)null)));
            Assert.That(quiet.Messages, Is.Empty);

            var normalizing = ExperimentalDesignCheck.Check(folder, spectra, normalizationRequested: true);
            Assert.That((normalizing.Verdict, normalizing.Problem), Is.EqualTo((ExperimentalDesignVerdict.Refuse, ExperimentalDesignProblem.MissingForNormalization)));
            Assert.That(normalizing.Messages.Single(), Does.Contain("Normalization requires a design (ExperimentalDesign.tsv or TmtDesign.txt)"));
        }

        /// <summary>The GUI refused this run ("Experimental design must be defined"): its only design is TmtDesign.txt.</summary>
        [Test]
        public static void ANormalizingRunWhoseDesignIsTmtDesignProceeds()
        {
            string folder = NewFolder("TmtOnly");
            string spectra = Path.Combine(folder, "file.1.mzML");
            WriteTmtDesign(folder, spectra, valid: true);

            var check = ExperimentalDesignCheck.Check(folder, new List<string> { spectra }, normalizationRequested: true);

            Assert.That(check.Verdict, Is.EqualTo(ExperimentalDesignVerdict.Proceed));
            Assert.That(check.DesignFileName, Is.EqualTo(GlobalVariables.TmtExperimentalDesignFileName));
            Assert.That(check.Messages.First(), Is.EqualTo("Read " + GlobalVariables.TmtExperimentalDesignFileName + " successfully"));
        }

        /// <summary>The GUI never read TmtDesign.txt, so it started a normalizing run on one it could not read.</summary>
        [Test]
        public static void AnUnreadableTmtDesignStopsANormalizingRunWithEveryError()
        {
            string folder = NewFolder("TmtBadNormalizing");
            string spectra = Path.Combine(folder, "file.1.mzML");
            WriteTmtDesign(folder, spectra, valid: false);

            var check = ExperimentalDesignCheck.Check(folder, new List<string> { spectra }, normalizationRequested: true);

            Assert.That((check.Verdict, check.Problem), Is.EqualTo((ExperimentalDesignVerdict.Refuse, ExperimentalDesignProblem.Unreadable)));
            Assert.That(check.Errors, Is.Not.Empty);
            Assert.That(check.Messages, Is.EqualTo(check.Errors));
            Assert.That(File.Exists(Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName)), "the file is never deleted");
        }

        [TestCase(true, "An experimental design file was found")]
        [TestCase(false, "A TMT design file was found")]
        public static void AnUnreadableDesignNothingNeedsAsksWhetherToContinue(bool classic, string questionStart)
        {
            string folder = NewFolder("BadQuiet" + (classic ? "Classic" : "Tmt"));
            string spectra = Path.Combine(folder, "file.1.mzML");
            if (classic) WriteClassicDesign(folder, spectra, valid: false);
            else WriteTmtDesign(folder, spectra, valid: false);

            var check = ExperimentalDesignCheck.Check(folder, new List<string> { spectra }, normalizationRequested: false);

            Assert.That((check.Verdict, check.Problem), Is.EqualTo((ExperimentalDesignVerdict.AskToContinueWithoutDesign, ExperimentalDesignProblem.Unreadable)));
            Assert.That(check.Question, Does.StartWith(questionStart).And.EndWith("Continue without an experimental design?"));
            Assert.That(check.Messages, Is.EqualTo(new[] { check.Question, "First error: " + check.Errors.First() }));
        }

        [Test]
        public static void AReadableDesignsWarningsFollowTheReadLine()
        {
            string folder = NewFolder("Warnings");
            string spectra = Path.Combine(folder, "file.mzML");
            WriteClassicDesign(folder, spectra, valid: true, 1, 3);
            var files = new List<string> { Path.Combine(folder, "file.mzML"), Path.Combine(folder, "file_1.mzML") };

            var check = ExperimentalDesignCheck.Check(folder, files, normalizationRequested: true);

            Assert.That(check.Verdict, Is.EqualTo(ExperimentalDesignVerdict.Proceed));
            Assert.That(check.Messages.First(), Is.EqualTo("Read " + GlobalVariables.ExperimentalDesignFileName + " successfully"));
            Assert.That(check.Messages.Skip(1).Single(), Does.Contain("biological replicate"));
        }
    }
}
