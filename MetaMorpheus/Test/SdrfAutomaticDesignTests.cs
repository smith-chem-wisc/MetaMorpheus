using EngineLayer;
using MetaMorpheusCommandLine;
using Nett;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// The automatic design (QuantProject gap AUTODESIGN, ruled 2026-11-07 and 2026-11-11): an SDRF beside the spectra files
    /// makes the run's experimental design and fills the search's own SDRF. These call
    /// <see cref="Program.ResolveSdrfDesign"/> directly for each branch, then drive <see cref="Program.Main"/> twice to prove
    /// the wiring: once refused, once through a whole search.
    /// </summary>
    [TestFixture]
    [NonParallelizable]
    internal static class SdrfAutomaticDesignTests
    {
        private const int Continue = 0;
        private const int Stop = 5;

        private static readonly string[] LabelFreeHeader =
        {
            "source name", "characteristics[organism]", "characteristics[organism part]", "characteristics[biological replicate]",
            "assay name", "comment[label]", "comment[fraction identifier]", "comment[technical replicate]", "comment[data file]",
            "factor value[genotype]",
        };

        private static string NewFolder(string name)
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "SdrfAutomaticDesign", name);
            if (Directory.Exists(folder)) Directory.Delete(folder, true);
            Directory.CreateDirectory(folder);
            return folder;
        }

        /// <summary>Empty spectra files in <paramref name="folder"/>, returned as full paths.</summary>
        private static List<string> Spectra(string folder, params string[] names)
        {
            foreach (string name in names)
            {
                File.WriteAllText(Path.Combine(folder, name), string.Empty);
            }
            return names.Select(n => Path.Combine(folder, n)).ToList();
        }

        /// <summary>
        /// A label-free SDRF beside the spectra: one row per (acquired file, genotype, biorep), naming the .raw file the
        /// search reads as .mzML, so the design is made through the join.
        /// </summary>
        private static string WriteSdrf(string folder, params (string RawFile, string Genotype, int Biorep)[] rows) =>
            WriteSdrf(folder, "deposit.sdrf.tsv", LabelFreeHeader, rows.Select(r => new[]
            {
                r.Genotype + " " + r.Biorep, "Mus musculus", "liver", r.Biorep.ToString(), "run " + r.RawFile, "label free sample",
                "1", "1", r.RawFile, r.Genotype,
            }).ToArray());

        private static string WriteSdrf(string folder, string name, string[] header, string[][] rows)
        {
            string path = Path.Combine(folder, name);
            File.WriteAllLines(path, new[] { string.Join("\t", header) }.Concat(rows.Select(r => string.Join("\t", r))));
            return path;
        }

        /// <summary>The same SDRF without its factor value column, so no condition can be read from it.</summary>
        private static string WriteSdrfWithNoFactorColumn(string folder, params (string RawFile, int Biorep)[] rows) =>
            WriteSdrf(folder, "deposit.sdrf.tsv", LabelFreeHeader.Take(LabelFreeHeader.Length - 1).ToArray(), rows.Select(r => new[]
            {
                "mouse " + r.Biorep, "Mus musculus", "liver", r.Biorep.ToString(), "run " + r.RawFile, "label free sample", "1", "1", r.RawFile,
            }).ToArray());

        private static List<(string, SearchTask)> SearchTasks(params SearchTask[] tasks) =>
            tasks.Select((t, i) => ("Task" + (i + 1) + "SearchTask", t)).ToList();

        private static int Resolve(string folder, List<string> spectra, List<(string, SearchTask)> tasks, List<string> console,
            bool normalize = false, bool reportToConsole = true, string answer = "n", CommandLineSettings settings = null) =>
            Resolve(folder, spectra, tasks, console, out _, normalize, reportToConsole, answer, settings);

        /// <summary>As above, and the SDRF for the run, which the command line hands to the runner.</summary>
        private static int Resolve(string folder, List<string> spectra, List<(string, SearchTask)> tasks, List<string> console,
            out SdrfForRun runSdrf, bool normalize = false, bool reportToConsole = true, string answer = "n",
            CommandLineSettings settings = null) =>
            Program.ResolveSdrfDesign(folder, spectra, tasks, settings ?? new CommandLineSettings(), normalize, reportToConsole,
                out runSdrf, console.Add, () => answer);

        /// <summary>The task's own SDRF settings, as the user saved them: the automatic design never changes them.</summary>
        private static void AssertTaskUnchanged(SearchTask task)
        {
            Assert.That(task.SearchParameters.InputSdrfPath, Is.Null, "the task's own InputSdrfPath is left as saved");
            Assert.That(task.SearchParameters.WriteSdrf, Is.False, "the task's own WriteSdrf is left as saved");
            Assert.That(task.RunSdrf, Is.Null, "only the runner hands a task the run's SDRF");
        }

        private static string DesignPath(string folder) => Path.Combine(folder, GlobalVariables.ExperimentalDesignFileName);

        // ---------------------------------------------------------------- nothing there, nothing changes

        [Test]
        public static void WithNoSdrfBesideTheSpectraNothingChanges()
        {
            string folder = NewFolder("NoSdrf");
            var spectra = Spectra(folder, "a1.mzML");
            var task = new SearchTask();
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(task), console, out var runSdrf), Is.EqualTo(Continue));
            Assert.That(Directory.GetFiles(folder).Select(Path.GetFileName), Is.EquivalentTo(new[] { "a1.mzML" }));
            AssertTaskUnchanged(task);
            Assert.That(runSdrf, Is.Null);
            Assert.That(console, Is.Empty);
        }

        // ---------------------------------------------------------------- the design and its marker

        /// <summary>
        /// The SDRF names X.raw and the search reads X.mzML: the design names the searched files, MetaMorpheus reads it back
        /// without error, a normalizing run's design check passes, the marker holds both sha256s, the SDRF comes back for the
        /// run with results.txt's line, the search's own SDRF is written from it, the task is not changed, and the SDRF is not
        /// touched.
        /// </summary>
        [Test]
        public static void AnSdrfBesideTheSpectraMakesTheDesignAndItsMarker()
        {
            string folder = NewFolder("Writes");
            var spectra = Spectra(folder, "a1.mzML", "b1.mzML");
            string sdrf = WriteSdrf(folder, ("a1.raw", "WT", 1), ("b1.raw", "KO", 1));
            string sdrfSha = SdrfDesignMarker.Sha256Of(sdrf);
            var task = new SearchTask();
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(task), console, out var runSdrf), Is.EqualTo(Continue), string.Join("\n", console));

            var readBack = ExperimentalDesign.ReadExperimentalDesign(DesignPath(folder), spectra, out var errors);
            Assert.That(errors, Is.Empty);
            Assert.That(readBack.Select(f => (f.FilenameWithoutExtension, f.Condition)), Is.EquivalentTo(new[] { ("a1", "WT"), ("b1", "KO") }));
            Assert.That(Program.ResolveExperimentalDesign(folder, spectra, normalizationRequested: true, reportToConsole: false), Is.EqualTo(Continue));

            var marker = SdrfDesignMarker.ReadFrom(folder);
            Assert.That(marker, Is.Not.Null);
            Assert.That(marker.DesignFileName, Is.EqualTo(GlobalVariables.ExperimentalDesignFileName));
            Assert.That(marker.DesignSha256, Is.EqualTo(SdrfDesignMarker.ContentSha256Of(DesignPath(folder))));
            Assert.That((marker.SdrfPath, marker.SdrfSha256), Is.EqualTo((Path.GetFullPath(sdrf), sdrfSha)));
            Assert.That(marker.Condition, Is.EqualTo("every factor value column, factor value[genotype]"));
            Assert.That(marker.StatesGrouping, Is.False, "nothing was stated, so the marker keeps no grouping");

            Assert.That(runSdrf, Is.EqualTo(new SdrfForRun(Path.GetFullPath(sdrf), marker.Describe())));
            AssertTaskUnchanged(task);
            Assert.That(SdrfForRun.InputSdrfPathFor(task.SearchParameters, runSdrf), Is.EqualTo(Path.GetFullPath(sdrf)));
            Assert.That(SdrfForRun.WritesSdrfFor(task.SearchParameters, runSdrf), Is.True);
            Assert.That(console, Has.Member("Task1SearchTask writes the search's own SDRF, because an SDRF is beside the spectra files."));
            Assert.That(SdrfDesignMarker.Sha256Of(sdrf), Is.EqualTo(sdrfSha), "the SDRF is never edited");
            Assert.That(console, Has.Some.StartWith("Wrote " + DesignPath(folder) + " from the SDRF"));
            Assert.That(string.Join("\n", console), Does.Not.Contain("plex"), "a label-free design has no plex to speak of");
        }

        /// <summary>
        /// A design the marker shows MetaMorpheus wrote, unedited, is written again for the run's files: here a later run
        /// searches one of the two files, and the design describes that one.
        /// </summary>
        [Test]
        public static void OurUneditedDesignIsWrittenAgainForTheRunsFiles()
        {
            string folder = NewFolder("Again");
            var spectra = Spectra(folder, "a1.mzML", "b1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1), ("b1.raw", "KO", 1));
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), new List<string>()), Is.EqualTo(Continue));

            var console = new List<string>();
            var onlyA1 = spectra.Take(1).ToList();
            Assert.That(Resolve(folder, onlyA1, SearchTasks(new SearchTask()), console), Is.EqualTo(Continue), string.Join("\n", console));

            var readBack = ExperimentalDesign.ReadExperimentalDesign(DesignPath(folder), onlyA1, out var errors);
            Assert.That(errors, Is.Empty);
            Assert.That(File.ReadAllLines(DesignPath(folder)), Has.Length.EqualTo(2), "header and a1 only");
            Assert.That(SdrfDesignMarker.ReadFrom(folder).DescribesUneditedDesignIn(folder), Is.True);
            Assert.That(console.Last(), Does.EndWith("again for this run's spectra files."));
        }

        /// <summary>QP-S41: an edited design is the user's, and is never touched; the SDRF is not read for it.</summary>
        [Test]
        public static void AnEditedDesignIsTheUsersAndIsNeverTouched()
        {
            string folder = NewFolder("Edited");
            var spectra = Spectra(folder, "a1.mzML", "b1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1), ("b1.raw", "KO", 1));
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), new List<string>()), Is.EqualTo(Continue));
            File.WriteAllText(DesignPath(folder), File.ReadAllText(DesignPath(folder)).Replace("KO", "knockout"));
            string edited = File.ReadAllText(DesignPath(folder));

            var task = new SearchTask();
            var console = new List<string>();
            Assert.That(Resolve(folder, spectra, SearchTasks(task), console, out var runSdrf), Is.EqualTo(Continue));

            Assert.That(File.ReadAllText(DesignPath(folder)), Is.EqualTo(edited));
            Assert.That(console, Has.Some.Contains("was not written by MetaMorpheus from an SDRF, or was edited since, so it is the design"));
            Assert.That(runSdrf.Provenance, Is.Null, "no design was written from the SDRF");
            Assert.That(runSdrf.SdrfPath, Is.Not.Null, "the SDRF still fills the search's own SDRF");
            AssertTaskUnchanged(task);
        }

        /// <summary>
        /// Ruled 2026-11-13: only a real edit makes a design the user's. The GUI's design window saves Windows line endings,
        /// so a design saved unchanged differs from ours byte for byte; it is still ours, and written again for the run.
        /// </summary>
        [Test]
        public static void OurDesignSavedWithOtherLineEndingsIsStillOurs()
        {
            string folder = NewFolder("Resaved");
            var spectra = Spectra(folder, "a1.mzML", "b1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1), ("b1.raw", "KO", 1));
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), new List<string>()), Is.EqualTo(Continue));
            string ourSha = SdrfDesignMarker.Sha256Of(DesignPath(folder));
            File.WriteAllLines(DesignPath(folder), File.ReadAllLines(DesignPath(folder)));
            Assume.That(SdrfDesignMarker.Sha256Of(DesignPath(folder)), Is.Not.EqualTo(ourSha), "the bytes changed");

            var console = new List<string>();
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), console), Is.EqualTo(Continue));

            Assert.That(console.Last(), Does.EndWith("again for this run's spectra files."));
            Assert.That(SdrfDesignMarker.ReadFrom(folder).DescribesUneditedDesignIn(folder), Is.True);
        }

        [Test]
        public static void ADesignWithNoMarkerIsTheUsersAndIsNeverTouched()
        {
            string folder = NewFolder("Users");
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1));
            File.WriteAllLines(DesignPath(folder), new[] { "FileName\tCondition\tBiorep\tFraction\tTechrep", "a1.mzML\tmine\t1\t1\t1" });
            string users = File.ReadAllText(DesignPath(folder));

            var console = new List<string>();
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), console), Is.EqualTo(Continue));

            Assert.That(File.ReadAllText(DesignPath(folder)), Is.EqualTo(users));
            Assert.That(File.Exists(SdrfDesignMarker.PathIn(folder)), Is.False);
        }

        // ---------------------------------------------------------------- a design that cannot be made

        /// <summary>QP-S42 (b): no factor column is refused; normalizing, the run stops with every reason and writes nothing.</summary>
        [Test]
        public static void ARefusedSdrfStopsANormalizingRunAndWritesNothing()
        {
            string folder = NewFolder("RefusedNormalizing");
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrfWithNoFactorColumn(folder, ("a1.raw", 1));
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), console, normalize: true), Is.EqualTo(Stop));

            Assert.That(File.Exists(DesignPath(folder)), Is.False);
            Assert.That(File.Exists(SdrfDesignMarker.PathIn(folder)), Is.False);
            Assert.That(string.Join("\n", console), Does.Contain("No condition")
                .And.Contain("On the command line, SingleCondition is --sdrfSingleCondition.")
                .And.Contain("Normalization needs an experimental design, so the run stops."));
        }

        /// <summary>As ResolveExperimentalDesign's prompt for an unreadable design: y or yes continues, anything else stops.</summary>
        [TestCase("y", Continue)]
        [TestCase("YES", Continue)]
        [TestCase("n", Stop)]
        [TestCase("", Stop)]
        [TestCase(null, Stop)]
        public static void ARefusedSdrfAsksWhetherToContinueWithoutADesign(string answer, int expected)
        {
            string folder = NewFolder("RefusedAsks");
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrfWithNoFactorColumn(folder, ("a1.raw", 1));
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), console, answer: answer), Is.EqualTo(expected));
            Assert.That(console.Last(), Does.EndWith("Continue without an experimental design? y/n"));
            Assert.That(File.Exists(DesignPath(folder)), Is.False);
        }

        [Test]
        public static void ARefusedSdrfAtVerbosityNoneContinuesWithoutADesign()
        {
            string folder = NewFolder("RefusedSilent");
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrfWithNoFactorColumn(folder, ("a1.raw", 1));
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), console, reportToConsole: false), Is.EqualTo(Continue));
            Assert.That(console, Is.Empty);
            Assert.That(File.Exists(DesignPath(folder)), Is.False);
        }

        /// <summary>
        /// A design MetaMorpheus wrote from the SDRF before, unedited, is removed with its marker when the SDRF no longer
        /// gives one: the run must not quantify with a design the SDRF no longer supports.
        /// </summary>
        [Test]
        public static void OurOldDesignIsRemovedWhenTheSdrfNoLongerGivesOne()
        {
            string folder = NewFolder("RemovedStale");
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1));
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), new List<string>()), Is.EqualTo(Continue));
            Assert.That(File.Exists(DesignPath(folder)));

            WriteSdrfWithNoFactorColumn(folder, ("a1.raw", 1));
            var console = new List<string>();
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), console, answer: "y"), Is.EqualTo(Continue));

            Assert.That(File.Exists(DesignPath(folder)), Is.False);
            Assert.That(File.Exists(SdrfDesignMarker.PathIn(folder)), Is.False);
            Assert.That(console, Has.Some.EndsWith("was removed: MetaMorpheus wrote it from the SDRF before, and it no longer fits."));
        }

        /// <summary>
        /// A design MetaMorpheus wrote that it cannot write again (read-only here, as on a data share) is refused, never
        /// thrown: a normalizing run stops with the reason. The design cannot be removed either, so its marker stays and the
        /// next run still knows it is MetaMorpheus's, rather than taking it for the user's.
        /// </summary>
        [Test]
        public static void ADesignThatCannotBeWrittenIsRefusedNotThrown()
        {
            string folder = NewFolder("ReadOnlyDesign");
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1));
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), new List<string>()), Is.EqualTo(Continue));
            File.SetAttributes(DesignPath(folder), FileAttributes.ReadOnly);
            try
            {
                var console = new List<string>();
                Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), console, normalize: true), Is.EqualTo(Stop));

                Assert.That(console, Has.Some.StartWith(DesignPath(folder) + " could not be removed:"));
                Assert.That(string.Join("\n", console), Does.Contain("but the experimental design made from it could not be written there:"));
                Assert.That(File.Exists(SdrfDesignMarker.PathIn(folder)), "the marker stays while the design does");
            }
            finally
            {
                File.SetAttributes(DesignPath(folder), FileAttributes.Normal);
            }
        }

        /// <summary>
        /// A design MetaMorpheus wrote that no longer fits and cannot be removed (read-only here) stops the run even when
        /// nothing normalizes, at verbosity none too, and asks nothing: a run that went on "without a design" would read it
        /// and quantify with the old conditions (Alexander-Sol's review of #2980). Both ways a design stops fitting: the SDRF
        /// now gives another design (WT became KO), which cannot be written over it, or gives none.
        /// </summary>
        [TestCase(true, true)]
        [TestCase(true, false)]
        [TestCase(false, true)]
        [TestCase(false, false)]
        public static void AnOldDesignThatCannotBeRemovedStopsTheRunWhateverItsNormalization(bool sdrfGivesAnotherDesign,
            bool reportToConsole)
        {
            string folder = NewFolder("StaleReadOnly_" + sdrfGivesAnotherDesign + reportToConsole);
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1));
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), new List<string>()), Is.EqualTo(Continue));
            string old = File.ReadAllText(DesignPath(folder));
            if (sdrfGivesAnotherDesign)
            {
                WriteSdrf(folder, ("a1.raw", "KO", 1));
            }
            else
            {
                WriteSdrfWithNoFactorColumn(folder, ("a1.raw", 1));
            }
            File.SetAttributes(DesignPath(folder), FileAttributes.ReadOnly);
            try
            {
                Assert.That(SdrfAutomaticDesign.Resolve(folder, spectra, SearchTasks(new SearchTask()), null, null).StopsTheRun, Is.True);

                var console = new List<string>();
                bool asked = false;
                Assert.That(Program.ResolveSdrfDesign(folder, spectra, SearchTasks(new SearchTask()), new CommandLineSettings(),
                        normalizationRequested: false, reportToConsole, out _, console.Add, () => { asked = true; return "y"; }),
                    Is.EqualTo(Stop), string.Join("\n", console));

                Assert.That(asked, Is.False, "nothing to ask: going on would quantify with the old design");
                Assert.That(File.ReadAllText(DesignPath(folder)), Is.EqualTo(old), "the old design is still there");
                Assert.That(File.Exists(SdrfDesignMarker.PathIn(folder)), "the marker stays while the design does");
                if (reportToConsole)
                {
                    Assert.That(console, Has.Some.StartWith(DesignPath(folder) + " could not be removed:"));
                    Assert.That(console.Last(), Does.EndWith("a run that went on would quantify with it, so the run stops."));
                }
                else
                {
                    Assert.That(console, Is.Empty);
                }
            }
            finally
            {
                File.SetAttributes(DesignPath(folder), FileAttributes.Normal);
            }
        }

        /// <summary>
        /// A marker that cannot be removed once its design is gone does not stop the run: nothing reads a marker without its
        /// design, so a run that goes on has no design, as the user agreed.
        /// </summary>
        [Test]
        public static void AMarkerLeftWithoutItsDesignDoesNotStopTheRun()
        {
            string folder = NewFolder("ReadOnlyMarker");
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1));
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), new List<string>()), Is.EqualTo(Continue));
            WriteSdrfWithNoFactorColumn(folder, ("a1.raw", 1));
            File.SetAttributes(SdrfDesignMarker.PathIn(folder), FileAttributes.ReadOnly);
            try
            {
                var console = new List<string>();
                Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), console, answer: "y"), Is.EqualTo(Continue));

                Assert.That(File.Exists(DesignPath(folder)), Is.False);
                Assert.That(console, Has.Some.StartWith(SdrfDesignMarker.PathIn(folder) + " could not be removed:"));
                Assert.That(console.Last(), Does.EndWith("Continue without an experimental design? y/n"));
            }
            finally
            {
                File.SetAttributes(SdrfDesignMarker.PathIn(folder), FileAttributes.Normal);
            }
        }

        /// <summary>
        /// A refused run removes the design MetaMorpheus wrote, and with its marker the grouping it kept; the run says so,
        /// rather than letting the next run quietly go back to every factor value column.
        /// </summary>
        [Test]
        public static void ARefusalSaysTheKeptGroupingWentWithTheMarker()
        {
            string folder = NewFolder("KeptGroupingGone");
            var spectra = Spectra(folder, "a1.mzML", "b1.mzML");
            WriteSdrfWithNoFactorColumn(folder, ("a1.raw", 1), ("b1.raw", 2));
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), new List<string>(),
                settings: new CommandLineSettings { SdrfSingleCondition = "survey" }), Is.EqualTo(Continue));

            var withAStrayFile = spectra.Concat(Spectra(folder, "c1.mzML")).ToList(); // c1 has no row in the SDRF
            var console = new List<string>();
            Assert.That(Resolve(folder, withAStrayFile, SearchTasks(new SearchTask()), console), Is.EqualTo(Stop));

            Assert.That(File.Exists(SdrfDesignMarker.PathIn(folder)), Is.False);
            Assert.That(console, Has.Some.EndsWith("was removed: MetaMorpheus wrote it from the SDRF before, and it no longer fits. "
                + "The grouping kept with it (one condition, 'survey') went with its marker; state it again for the next run."));
        }

        /// <summary>Which SDRF, when a folder holds two, is not ruled, so none is chosen and none is handed on.</summary>
        [Test]
        public static void TwoSdrfsAreNotChosenBetween()
        {
            string folder = NewFolder("Two");
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1));
            File.Copy(Path.Combine(folder, "deposit.sdrf.tsv"), Path.Combine(folder, "other.sdrf.tsv"));
            var task = new SearchTask();
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(task), console, out var runSdrf, normalize: true), Is.EqualTo(Stop));
            Assert.That(console.Single(), Does.StartWith("2 SDRF files are beside the spectra files (deposit.sdrf.tsv, other.sdrf.tsv)"));
            Assert.That(runSdrf, Is.Null);
            AssertTaskUnchanged(task);
            Assert.That(File.Exists(DesignPath(folder)), Is.False);
        }

        /// <summary>
        /// Ruled 2026-11-13: a name ending .sdrf.tsv is an SDRF, and so is a file named just sdrf.tsv, case ignored. A name
        /// that only contains "sdrf" is not, so nothing is made from it and nothing is said.
        /// </summary>
        [TestCase("deposit.sdrf.tsv", true)]
        [TestCase("sdrf.tsv", true)]
        [TestCase("SDRF.TSV", true)]
        [TestCase("my_sdrf.tsv", false)]
        [TestCase("sdrf.tsv.bak", false)]
        [TestCase("._deposit.sdrf.tsv", false)]
        public static void AnSdrfIsKnownByItsFileName(string name, bool isSdrf)
        {
            // Windows folder names ignore case, so the upper-case name gets a folder of its own.
            string folder = NewFolder("Named_" + name.Replace('.', '_') + (name.Any(char.IsUpper) ? "_upper" : ""));
            var spectra = Spectra(folder, "a1.mzML");
            string written = WriteSdrf(folder, ("a1.raw", "WT", 1));
            string renamed = Path.Combine(folder, name);
            File.Move(written, renamed);
            var task = new SearchTask();
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(task), console, out var runSdrf), Is.EqualTo(Continue), string.Join("\n", console));

            Assert.That(File.Exists(DesignPath(folder)), Is.EqualTo(isSdrf));
            Assert.That(runSdrf?.SdrfPath, Is.EqualTo(isSdrf ? Path.GetFullPath(renamed) : null));
            if (!isSdrf)
            {
                Assert.That(console, Is.Empty);
            }
        }

        /// <summary>A file named just sdrf.tsv beside a named SDRF is two SDRFs: neither is chosen.</summary>
        [Test]
        public static void ABareSdrfBesideANamedOneIsNotChosenBetween()
        {
            string folder = NewFolder("TwoWithBare");
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1));
            File.Copy(Path.Combine(folder, "deposit.sdrf.tsv"), Path.Combine(folder, "sdrf.tsv"));
            var task = new SearchTask();
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(task), console, out var runSdrf, normalize: true), Is.EqualTo(Stop));
            Assert.That(console.Single(), Does.StartWith("2 SDRF files are beside the spectra files (deposit.sdrf.tsv, sdrf.tsv)"));
            Assert.That(runSdrf, Is.Null);
            Assert.That(File.Exists(DesignPath(folder)), Is.False);
        }

        /// <summary>The AppleDouble file macOS leaves beside a copied SDRF is not a second SDRF: the design is made from the real one.</summary>
        [Test]
        public static void AMacOsCopyBesideTheSdrfIsNotASecondSdrf()
        {
            string folder = NewFolder("AppleDouble");
            var spectra = Spectra(folder, "a1.mzML");
            string sdrf = WriteSdrf(folder, ("a1.raw", "WT", 1));
            File.WriteAllBytes(Path.Combine(folder, "._deposit.sdrf.tsv"), new byte[] { 0, 5, 22, 7, 0, 2, 0, 0, 77, 97, 99, 32, 79, 83, 32, 88 });
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), console, out var runSdrf, normalize: true), Is.EqualTo(Continue),
                string.Join("\n", console));
            Assert.That(File.Exists(DesignPath(folder)));
            Assert.That(runSdrf.SdrfPath, Is.EqualTo(Path.GetFullPath(sdrf)));
        }

        /// <summary>
        /// Two SDRFs are refused only where a design would be made from one. Beside the user's own design, that design is read
        /// as with no SDRF, and neither SDRF goes to the search.
        /// </summary>
        [Test]
        public static void TwoSdrfsBesideTheUsersDesignLeaveItTheDesign()
        {
            string folder = NewFolder("TwoWithUsers");
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1));
            File.Copy(Path.Combine(folder, "deposit.sdrf.tsv"), Path.Combine(folder, "other.sdrf.tsv"));
            File.WriteAllLines(DesignPath(folder), new[] { "FileName\tCondition\tBiorep\tFraction\tTechrep", "a1.mzML\tmine\t1\t1\t1" });
            string users = File.ReadAllText(DesignPath(folder));
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), console, out var runSdrf, normalize: true), Is.EqualTo(Continue),
                string.Join("\n", console));
            Assert.That(File.ReadAllText(DesignPath(folder)), Is.EqualTo(users));
            Assert.That(console.Single(), Does.EndWith("so it is the design; the 2 SDRF files beside it are not read."));
            Assert.That(runSdrf, Is.Null);
        }

        /// <summary>Two SDRFs beside a design MetaMorpheus wrote: it no longer has one SDRF behind it, so it is removed with its marker.</summary>
        [Test]
        public static void TwoSdrfsRemoveTheDesignMetaMorpheusWrote()
        {
            string folder = NewFolder("TwoWithOurs");
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1));
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), new List<string>()), Is.EqualTo(Continue));
            File.Copy(Path.Combine(folder, "deposit.sdrf.tsv"), Path.Combine(folder, "other.sdrf.tsv"));
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), console, out var runSdrf, answer: "y"), Is.EqualTo(Continue));
            Assert.That(File.Exists(DesignPath(folder)), Is.False);
            Assert.That(File.Exists(SdrfDesignMarker.PathIn(folder)), Is.False);
            Assert.That(console, Has.Some.EndsWith("was removed: MetaMorpheus wrote it from the SDRF before, and it no longer fits."));
            Assert.That(console.Last(), Does.StartWith("2 SDRF files are beside the spectra files (deposit.sdrf.tsv, other.sdrf.tsv)"));
            Assert.That(runSdrf, Is.Null);
        }

        /// <summary>A search that quantifies nothing reads no design, so two SDRFs beside it stop nothing; neither goes to the search.</summary>
        [Test]
        public static void TwoSdrfsBesideASearchThatQuantifiesNothingStopNothing()
        {
            string folder = NewFolder("TwoNoQuant");
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1));
            File.Copy(Path.Combine(folder, "deposit.sdrf.tsv"), Path.Combine(folder, "other.sdrf.tsv"));
            var task = new SearchTask();
            task.SearchParameters.DoLabelFreeQuantification = false;
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(task), console, out var runSdrf, normalize: true), Is.EqualTo(Continue));
            Assert.That(console.Single(), Is.EqualTo("2 SDRF files are beside the spectra files (deposit.sdrf.tsv, other.sdrf.tsv), "
                + "and MetaMorpheus does not choose between them, so neither is given to the search."));
            Assert.That(runSdrf, Is.Null);
        }

        // ---------------------------------------------------------------- the grouping the user states

        private static List<string> Conditions(string folder, List<string> spectra) =>
            ExperimentalDesign.ReadExperimentalDesign(DesignPath(folder), spectra, out _).Select(f => f.Condition).ToList();

        private const string UsedAgain = "The grouping stated when MetaMorpheus last wrote ExperimentalDesign.tsv from the SDRF is used again";

        /// <summary>
        /// QP-S42 (b): the command line states one group; it is recorded as stated, and the marker keeps it as a key with the
        /// design it made. The task is not changed.
        /// </summary>
        [Test]
        public static void TheCommandLineStatesOneConditionAndTheMarkerKeepsIt()
        {
            string folder = NewFolder("StatedSingle");
            var spectra = Spectra(folder, "a1.mzML", "b1.mzML");
            WriteSdrfWithNoFactorColumn(folder, ("a1.raw", 1), ("b1.raw", 2));
            var task = new SearchTask();
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(task), console, out var runSdrf,
                settings: new CommandLineSettings { SdrfSingleCondition = "survey" }), Is.EqualTo(Continue), string.Join("\n", console));

            Assert.That(Conditions(folder, spectra), Is.All.EqualTo("survey"));
            var marker = SdrfDesignMarker.ReadFrom(folder);
            Assert.That(marker.Condition, Is.EqualTo("one condition, 'survey', stated by the user"));
            Assert.That((marker.SingleCondition, marker.ConditionColumns.Count, marker.UnknownFactorIsOwnCondition), Is.EqualTo(("survey", 0, false)));
            Assert.That(runSdrf.Provenance, Does.Contain("stated by the user"));
            Assert.That(console, Has.None.StartWith(UsedAgain), "nothing was kept from before");
            AssertTaskUnchanged(task);
        }

        /// <summary>
        /// The grouping a run stated is used again by a later run that states none, and the run says so: here the SDRF has no
        /// factor column, so without the kept grouping the second run would be refused.
        /// </summary>
        [Test]
        public static void AGroupingStatedBeforeIsUsedAgainAndSaidSo()
        {
            string folder = NewFolder("UsedAgain");
            var spectra = Spectra(folder, "a1.mzML", "b1.mzML");
            WriteSdrfWithNoFactorColumn(folder, ("a1.raw", 1), ("b1.raw", 2));
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), new List<string>(),
                settings: new CommandLineSettings { SdrfSingleCondition = "survey" }), Is.EqualTo(Continue));

            var console = new List<string>();
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), console), Is.EqualTo(Continue), string.Join("\n", console));

            Assert.That(Conditions(folder, spectra), Is.All.EqualTo("survey"));
            Assert.That(console, Has.Member(UsedAgain + " (one condition, 'survey'), as DesignFromSdrf.txt records it. Give --sdrfCondition, "
                + "--sdrfSingleCondition or --sdrfUnknownFactorIsCondition to state another, or delete the design file to make it again "
                + "from every factor value column."));
            Assert.That(SdrfDesignMarker.ReadFrom(folder).SingleCondition, Is.EqualTo("survey"), "kept for the run after");
        }

        /// <summary>
        /// The kept grouping is described as it was stated: columns, and files whose factor nobody stated as their own
        /// condition. Without it the second run is refused, because b1's factor is 'not available'.
        /// </summary>
        [TestCase(new[] { "genotype" }, "columns genotype, unknown factors their own condition")]
        [TestCase(new string[0], "every factor value column, unknown factors their own condition")]
        public static void AKeptGroupingIsDescribedAsStated(string[] columns, string described)
        {
            string folder = NewFolder("UsedAgainDescribed_" + columns.Length);
            var spectra = Spectra(folder, "a1.mzML", "b1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1), ("b1.raw", "not available", 1));
            var first = new List<string>();
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), first,
                settings: new CommandLineSettings { SdrfConditionColumns = columns, SdrfUnknownFactorIsCondition = true }),
                Is.EqualTo(Continue), string.Join("\n", first));

            var console = new List<string>();
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), console), Is.EqualTo(Continue), string.Join("\n", console));

            Assert.That(Conditions(folder, spectra), Is.EquivalentTo(new[] { "WT", "not available" }));
            Assert.That(console, Has.Some.StartWith(UsedAgain + " (" + described + ")"));
        }

        /// <summary>
        /// Two samples, both biological replicate 1, with two factor columns: genotype (WT, KO) and treatment (DMSO for
        /// both). Every grouping these tests state keeps the two samples in different conditions; one condition, or
        /// treatment alone, would make them duplicates, which is refused (QP-S45).
        /// </summary>
        private static void WriteTwoFactorSdrf(string folder) =>
            WriteSdrf(folder, "deposit.sdrf.tsv", LabelFreeHeader.Append("factor value[treatment]").ToArray(), new[]
            {
                new[] { "WT 1", "Mus musculus", "liver", "1", "run a1.raw", "label free sample", "1", "1", "a1.raw", "WT", "DMSO" },
                new[] { "KO 1", "Mus musculus", "liver", "1", "run b1.raw", "label free sample", "1", "1", "b1.raw", "KO", "DMSO" },
            });

        /// <summary>A grouping stated now replaces the one kept, and the marker keeps the new one.</summary>
        [Test]
        public static void AStatedGroupingReplacesTheOneKept()
        {
            string folder = NewFolder("StatedReplaces");
            var spectra = Spectra(folder, "a1.mzML", "b1.mzML");
            WriteTwoFactorSdrf(folder);
            var first = new List<string>();
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), first,
                settings: new CommandLineSettings { SdrfConditionColumns = new[] { "treatment", "genotype" } }), Is.EqualTo(Continue),
                string.Join("\n", first));
            Assert.That(Conditions(folder, spectra), Has.All.Contains("DMSO"), "both columns make the condition");

            var console = new List<string>();
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), console,
                settings: new CommandLineSettings { SdrfConditionColumns = new[] { "genotype" } }), Is.EqualTo(Continue), string.Join("\n", console));

            Assert.That(Conditions(folder, spectra), Is.EquivalentTo(new[] { "WT", "KO" }));
            Assert.That(SdrfDesignMarker.ReadFrom(folder).ConditionColumns, Is.EqualTo(new[] { "genotype" }));
            Assert.That(console, Has.None.StartWith(UsedAgain));
        }

        /// <summary>
        /// Deleting the design is the way back to every factor value column: a marker beside no design keeps nothing for the
        /// next run, which makes the design afresh.
        /// </summary>
        [Test]
        public static void DeletingTheDesignMakesItAgainFromEveryFactorColumn()
        {
            string folder = NewFolder("DeletedDesign");
            var spectra = Spectra(folder, "a1.mzML", "b1.mzML");
            WriteTwoFactorSdrf(folder);
            var first = new List<string>();
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), first,
                settings: new CommandLineSettings { SdrfConditionColumns = new[] { "genotype" } }), Is.EqualTo(Continue), string.Join("\n", first));
            Assert.That(Conditions(folder, spectra), Is.EquivalentTo(new[] { "WT", "KO" }));
            File.Delete(DesignPath(folder));
            Assume.That(SdrfDesignMarker.ReadFrom(folder).ConditionColumns, Is.EqualTo(new[] { "genotype" }), "the marker is still there");

            var console = new List<string>();
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), console), Is.EqualTo(Continue), string.Join("\n", console));

            Assert.That(Conditions(folder, spectra), Has.All.Contains("DMSO"), "every factor value column, treatment included");
            Assert.That(SdrfDesignMarker.ReadFrom(folder).StatesGrouping, Is.False);
            Assert.That(console, Has.None.StartWith(UsedAgain));
        }

        // ---------------------------------------------------------------- the search's own SDRF, and what reads no design

        [Test]
        public static void ATasksOwnInputSdrfIsKeptAndSaidToDiffer()
        {
            string folder = NewFolder("OwnInput");
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1));
            string own = Path.Combine(NewFolder("OwnInputElsewhere"), "curated.sdrf.tsv");
            var task = new SearchTask();
            task.SearchParameters.InputSdrfPath = own;
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(task), console, out var runSdrf), Is.EqualTo(Continue));
            Assert.That(task.SearchParameters.InputSdrfPath, Is.EqualTo(own));
            Assert.That(SdrfForRun.InputSdrfPathFor(task.SearchParameters, runSdrf), Is.EqualTo(own), "the task's own SDRF wins");
            Assert.That(console, Has.Some.StartWith($"Task1SearchTask names its own InputSdrfPath ({own}), not the SDRF beside the spectra files"));
        }

        [Test]
        public static void ASearchThatQuantifiesNothingGetsNoDesignButStillItsSdrf()
        {
            string folder = NewFolder("NoQuant");
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrfWithNoFactorColumn(folder, ("a1.raw", 1));
            var task = new SearchTask();
            task.SearchParameters.DoLabelFreeQuantification = false;
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(task), console, out var runSdrf, normalize: true), Is.EqualTo(Continue));
            Assert.That(File.Exists(DesignPath(folder)), Is.False);
            Assert.That(runSdrf.SdrfPath, Is.Not.Null);
            Assert.That(runSdrf.Provenance, Is.Null);
            Assert.That(SdrfForRun.WritesSdrfFor(task.SearchParameters, runSdrf), Is.True);
        }

        /// <summary>
        /// The values a search uses for the run: the task's own input SDRF wins over the run's, and the search writes its
        /// own SDRF when the task asks for it or the run has an SDRF. Without parameters or a run SDRF, neither.
        /// </summary>
        [Test]
        public static void TheRunsSdrfFillsOnlyWhatTheTaskLeftBlank()
        {
            var run = new SdrfForRun(@"C:\deposit\deposit.sdrf.tsv", null);
            var blank = new SearchParameters { InputSdrfPath = "  " };
            var own = new SearchParameters { InputSdrfPath = @"C:\curated.sdrf.tsv", WriteSdrf = true };

            Assert.That(SdrfForRun.InputSdrfPathFor(null, null), Is.Null);
            Assert.That(SdrfForRun.InputSdrfPathFor(blank, null), Is.Null);
            Assert.That(SdrfForRun.InputSdrfPathFor(blank, run), Is.EqualTo(run.SdrfPath));
            Assert.That(SdrfForRun.InputSdrfPathFor(own, run), Is.EqualTo(own.InputSdrfPath));

            Assert.That(SdrfForRun.WritesSdrfFor(null, null), Is.False);
            Assert.That(SdrfForRun.WritesSdrfFor(blank, null), Is.False);
            Assert.That(SdrfForRun.WritesSdrfFor(blank, run), Is.True);
            Assert.That(SdrfForRun.WritesSdrfFor(own, null), Is.True);
        }

        // ---------------------------------------------------------------- isobaric

        private static string WriteTmtSdrf(string folder) => WriteSdrf(folder, "deposit.sdrf.tsv",
            new[] { "source name", "characteristics[biological replicate]", "comment[label]", "comment[data file]", "comment[plex]", "factor value[genotype]" },
            new[]
            {
                new[] { "S1", "1", "TMT126", "p1.raw", "P1", "WT" },
                new[] { "S2", "1", "TMT127N", "p1.raw", "P1", "KO" },
            });

        private static SearchTask Tmt11Search() => MultiplexSearch("TMT11");

        private static SearchTask MultiplexSearch(string label)
        {
            var task = new SearchTask();
            task.SearchParameters.DoMultiplexQuantification = true;
            task.SearchParameters.MultiplexModId = label;
            return task;
        }

        /// <summary>The kit is the search's multiplex label; the plex is read only from a stated source (QP-S43).</summary>
        [Test]
        public static void AnIsobaricSearchWritesTmtDesignFromAStatedPlex()
        {
            string folder = NewFolder("Tmt");
            var spectra = Spectra(folder, "p1.mzML");
            WriteTmtSdrf(folder);
            var console = new List<string>();

            int code = Resolve(folder, spectra, SearchTasks(Tmt11Search()), console, settings: new CommandLineSettings { SdrfPlexColumn = "comment[plex]" });

            Assert.That(code, Is.EqualTo(Continue), string.Join("\n", console));
            string tmtDesign = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);
            TmtExperimentalDesign.Read(tmtDesign, spectra, out var errors);
            Assert.That(errors, Is.Empty);
            Assert.That(SdrfDesignMarker.ReadFrom(folder).DesignFileName, Is.EqualTo(GlobalVariables.TmtExperimentalDesignFileName));
            Assert.That(File.Exists(DesignPath(folder)), Is.False);
        }

        /// <summary>Ruled 2026-11-13: with no plex stated, the SDRF's comment[plex] column is read.</summary>
        [Test]
        public static void AnIsobaricSdrfsCommentPlexIsReadWhenNoPlexIsStated()
        {
            string folder = NewFolder("TmtDefaultPlex");
            var spectra = Spectra(folder, "p1.mzML");
            WriteTmtSdrf(folder);
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(Tmt11Search()), console, normalize: true), Is.EqualTo(Continue),
                string.Join("\n", console));

            var files = TmtExperimentalDesign.Read(Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName), spectra, out var errors);
            Assert.That(errors, Is.Empty);
            Assert.That(files.Select(f => f.Plex).Distinct(), Is.EqualTo(new[] { "P1" }));
            // mzLib's report says the column was "declared"; the run says MetaMorpheus chose it.
            Assert.That(string.Join("\n", console),
                Does.Contain("No plex source was named, so MetaMorpheus read the plex from the SDRF's comment[plex] column."));
        }

        /// <summary>A stated plex column wins over comment[plex].</summary>
        [Test]
        public static void AStatedPlexColumnWinsOverCommentPlex()
        {
            string folder = NewFolder("TmtStatedPlex");
            var spectra = Spectra(folder, "p1.mzML");
            WriteSdrf(folder, "deposit.sdrf.tsv",
                new[] { "source name", "characteristics[biological replicate]", "comment[label]", "comment[data file]", "comment[plex]", "comment[batch]", "factor value[genotype]" },
                new[]
                {
                    new[] { "S1", "1", "TMT126", "p1.raw", "P1", "B9", "WT" },
                    new[] { "S2", "1", "TMT127N", "p1.raw", "P1", "B9", "KO" },
                });
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(Tmt11Search()), console,
                settings: new CommandLineSettings { SdrfPlexColumn = "comment[batch]" }), Is.EqualTo(Continue), string.Join("\n", console));

            var files = TmtExperimentalDesign.Read(Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName), spectra, out _);
            Assert.That(files.Select(f => f.Plex).Distinct(), Is.EqualTo(new[] { "B9" }));
            Assert.That(string.Join("\n", console), Does.Not.Contain("No plex source was named"));
        }

        /// <summary>No comment[plex] and no plex stated: refused, and the report says comment[plex] was looked for.</summary>
        [Test]
        public static void AnIsobaricSdrfWithNoCommentPlexAndNoPlexStatedIsRefused()
        {
            string folder = NewFolder("TmtNoPlex");
            var spectra = Spectra(folder, "p1.mzML");
            WriteSdrf(folder, "deposit.sdrf.tsv",
                new[] { "source name", "characteristics[biological replicate]", "comment[label]", "comment[data file]", "factor value[genotype]" },
                new[]
                {
                    new[] { "S1", "1", "TMT126", "p1.raw", "WT" },
                    new[] { "S2", "1", "TMT127N", "p1.raw", "KO" },
                });
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(Tmt11Search()), console, normalize: true), Is.EqualTo(Stop));

            Assert.That(string.Join("\n", console), Does.Contain("No plex source was declared")
                .And.Contain("The SDRF has no comment[plex] column, which MetaMorpheus reads for the plex when none is named."));
            Assert.That(File.Exists(Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName)), Is.False);
        }

        /// <summary>TMT16 is searched with MetaMorpheus's TMT18 label; a design made for a TMT16 search says so.</summary>
        [Test]
        public static void ATmt16DesignSaysItIsReadByATmt18Search()
        {
            string folder = NewFolder("Tmt16");
            var spectra = Spectra(folder, "p1.mzML");
            WriteTmtSdrf(folder);
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(MultiplexSearch("TMT16")), console), Is.EqualTo(Continue), string.Join("\n", console));

            Assert.That(File.Exists(Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName)));
            Assert.That(console, Has.Member(SdrfAutomaticDesign.Tmt16Note));
        }

        /// <summary>
        /// One design serves every search task, so a label-free task beside a TMT task is refused, each task named with how
        /// it quantifies.
        /// </summary>
        [Test]
        public static void SearchTasksThatQuantifyDifferentlyAreRefused()
        {
            string folder = NewFolder("TasksQuantifyDifferently");
            var spectra = Spectra(folder, "p1.mzML");
            WriteTmtSdrf(folder);
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask(), Tmt11Search()), console, normalize: true), Is.EqualTo(Stop));

            Assert.That(console.Last(), Does.StartWith("The search tasks quantify differently (Task1SearchTask: label-free, Task2SearchTask: TMT11), "
                + "and one experimental design cannot serve them all, so none was made from the SDRF."));
            Assert.That(File.Exists(DesignPath(folder)), Is.False);
            Assert.That(File.Exists(Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName)), Is.False);
        }

        /// <summary>
        /// A search task that quantifies nothing reads no design, so it is not counted as label-free beside a TMT task: the
        /// TMT design is made for the task that quantifies.
        /// </summary>
        [Test]
        public static void ASearchTaskThatQuantifiesNothingDoesNotBlockTheTmtDesign()
        {
            string folder = NewFolder("IdOnlyBesideTmt");
            var spectra = Spectra(folder, "p1.mzML");
            WriteTmtSdrf(folder);
            var identificationOnly = new SearchTask();
            identificationOnly.SearchParameters.DoLabelFreeQuantification = false;
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(identificationOnly, Tmt11Search()), console, normalize: true), Is.EqualTo(Continue),
                string.Join("\n", console));
            Assert.That(File.Exists(Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName)));
            Assert.That(string.Join("\n", console), Does.Not.Contain("quantify differently"));
        }

        /// <summary>A multiplex label mzLib has no isobaric kit for, or none at all, cannot read the SDRF's channels: refused.</summary>
        [TestCase("NoSuchKit")]
        [TestCase((string)null)]
        public static void AMultiplexLabelMzLibDoesNotKnowIsRefused(string label)
        {
            string folder = NewFolder("UnknownKit_" + (label ?? "none"));
            var spectra = Spectra(folder, "p1.mzML");
            WriteTmtSdrf(folder);
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(MultiplexSearch(label)), console, normalize: true), Is.EqualTo(Stop));

            Assert.That(console.Last(), Does.StartWith($"The search's multiplex label '{label}' is not an isobaric kit mzLib knows, "
                + "so the SDRF's channels cannot be read."));
            Assert.That(File.Exists(Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName)), Is.False);
        }

        /// <summary>
        /// Our unedited label-free design, when the search becomes a TMT search: it is removed and the TMT design written in
        /// its place, so the folder never holds both kinds (which stops a run).
        /// </summary>
        [Test]
        public static void OurDesignOfTheOtherKindIsReplacedWhenTheSearchsLabelChanges()
        {
            string folder = NewFolder("LabelChanged");
            var spectra = Spectra(folder, "p1.mzML");
            WriteSdrf(folder, ("p1.raw", "WT", 1));
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), new List<string>()), Is.EqualTo(Continue));
            Assert.That(File.Exists(DesignPath(folder)));

            WriteTmtSdrf(folder);
            var console = new List<string>();
            Assert.That(Resolve(folder, spectra, SearchTasks(Tmt11Search()), console), Is.EqualTo(Continue), string.Join("\n", console));

            string tmtDesign = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);
            Assert.That(File.Exists(DesignPath(folder)), Is.False);
            Assert.That(File.Exists(tmtDesign));
            Assert.That(SdrfDesignMarker.ReadFrom(folder).DesignFileName, Is.EqualTo(GlobalVariables.TmtExperimentalDesignFileName));
            Assert.That(console.Last(), Does.StartWith("Wrote " + tmtDesign + " from the SDRF").And.EndWith("again for this run's spectra files."));
        }

        // ---------------------------------------------------------------- Resolve itself, and the marker

        /// <summary>
        /// Resolve with nothing stated and the default wording (null for both): the outcome names the SDRF beside the spectra
        /// and the design it wrote.
        /// </summary>
        [Test]
        public static void TheOutcomeNamesTheSdrfAndTheDesignWritten()
        {
            string folder = NewFolder("Outcome");
            var spectra = Spectra(folder, "a1.mzML");
            string sdrf = WriteSdrf(folder, ("a1.raw", "WT", 1));

            var outcome = SdrfAutomaticDesign.Resolve(folder, spectra, SearchTasks(new SearchTask()), null, null);

            Assert.That(outcome.Refusal, Is.Null, string.Join("\n", outcome.Messages));
            Assert.That(outcome.SdrfPath, Is.EqualTo(Path.GetFullPath(sdrf)));
            Assert.That(outcome.DesignPath, Is.EqualTo(DesignPath(folder)));
            Assert.That(outcome.RunSdrf, Is.EqualTo(new SdrfForRun(Path.GetFullPath(sdrf), SdrfDesignMarker.ReadFrom(folder).Describe())));
        }

        /// <summary>
        /// With the default wording (no front end names its controls), a kept grouping is used again and the way back to
        /// every factor value column is named.
        /// </summary>
        [Test]
        public static void TheDefaultWordingSaysHowToGoBackToEveryFactorColumn()
        {
            string folder = NewFolder("DefaultWording");
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1));
            var stated = new SdrfDesignStatements { SingleCondition = "all" };
            Assume.That(SdrfAutomaticDesign.Resolve(folder, spectra, SearchTasks(new SearchTask()), stated, null).Refusal, Is.Null);

            var outcome = SdrfAutomaticDesign.Resolve(folder, spectra, SearchTasks(new SearchTask()), null, null);

            Assert.That(outcome.Messages, Has.Member(UsedAgain + " (one condition, 'all'), as DesignFromSdrf.txt records it. "
                + "To make it again from every factor value column, delete the design file."));
        }

        /// <summary>
        /// The marker keeps a stated grouping as keys, columns tab-separated, and reads it back; a marker stating nothing
        /// writes no grouping key, and a value other than 'true' states no unknown-factor condition.
        /// </summary>
        [Test]
        public static void TheMarkerKeepsTheStatedGroupingAsKeys()
        {
            string folder = NewFolder("MarkerGrouping");
            new SdrfDesignMarker(GlobalVariables.ExperimentalDesignFileName, "abc", "x.sdrf.tsv", "def", "described",
                new[] { "genotype", "factor value[treatment]" }, "one", true).WriteTo(folder);
            string text = File.ReadAllText(SdrfDesignMarker.PathIn(folder));

            Assert.That(text, Does.Contain("ConditionColumns\tgenotype\tfactor value[treatment]\n")
                .And.Contain("SingleCondition\tone\n").And.Contain("UnknownFactorIsOwnCondition\ttrue\n"));
            var read = SdrfDesignMarker.ReadFrom(folder);
            Assert.That(read.ConditionColumns, Is.EqualTo(new[] { "genotype", "factor value[treatment]" }));
            Assert.That((read.SingleCondition, read.UnknownFactorIsOwnCondition, read.StatesGrouping), Is.EqualTo(("one", true, true)));

            File.WriteAllText(SdrfDesignMarker.PathIn(folder), text.Replace("UnknownFactorIsOwnCondition\ttrue", "UnknownFactorIsOwnCondition\tyes"));
            Assert.That(SdrfDesignMarker.ReadFrom(folder).UnknownFactorIsOwnCondition, Is.False);

            new SdrfDesignMarker(GlobalVariables.ExperimentalDesignFileName, "abc", "x.sdrf.tsv", "def", "described").WriteTo(folder);
            string none = File.ReadAllText(SdrfDesignMarker.PathIn(folder));
            // The comment line names the keys, so look for each key followed by its tab.
            Assert.That(none, Does.Not.Contain("ConditionColumns\t").And.Not.Contain("SingleCondition\t").And.Not.Contain("UnknownFactorIsOwnCondition\t"));
            var nothing = SdrfDesignMarker.ReadFrom(folder);
            Assert.That((nothing.ConditionColumns.Count, nothing.SingleCondition, nothing.UnknownFactorIsOwnCondition, nothing.StatesGrouping),
                Is.EqualTo((0, (string)null, false, false)));
        }

        // ---------------------------------------------------------------- the runner

        /// <summary>
        /// The runner copies the marker of a design written from the SDRF into the run's output folder, beside Task Settings,
        /// from the folder of the spectra the run started from. With no design written from the SDRF (no provenance) or no
        /// SDRF for the run, nothing is copied; a marker that cannot be copied is warned about. No task runs here: the copy is
        /// the runner's own, after its tasks.
        /// </summary>
        [TestCase(true, true, true)]
        [TestCase(true, false, false)]
        [TestCase(false, true, false)]
        public static void TheRunnerCopiesTheMarkerIntoTheOutputFolder(bool withRunSdrf, bool withProvenance, bool copied)
        {
            string folder = NewFolder("RunnerCopies_" + withRunSdrf + withProvenance);
            var spectra = Spectra(folder, "a1.mzML");
            new SdrfDesignMarker(GlobalVariables.ExperimentalDesignFileName, "abc", "x.sdrf.tsv", "def", "described").WriteTo(folder);
            string output = Path.Combine(folder, "output");
            Directory.CreateDirectory(output);
            var runSdrf = withRunSdrf ? new SdrfForRun(Path.Combine(folder, "x.sdrf.tsv"), withProvenance ? "provenance" : null) : null;

            var runner = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)>(), spectra, new List<EngineLayer.DatabaseLoading.DbForTask>(),
                output, runSdrf);
            runner.Run();

            string copy = Path.Combine(output, SdrfDesignMarker.FileName);
            Assert.That(File.Exists(copy), Is.EqualTo(copied));
            if (copied)
            {
                Assert.That(File.ReadAllText(copy), Is.EqualTo(File.ReadAllText(SdrfDesignMarker.PathIn(folder))));
            }
            Assert.That(runner.Warnings, Is.Empty);
        }

        [Test]
        public static void AMarkerTheRunnerCannotCopyIsWarnedAbout()
        {
            string folder = NewFolder("RunnerCannotCopy");
            var spectra = Spectra(folder, "a1.mzML");
            string output = Path.Combine(folder, "output");
            Directory.CreateDirectory(output);

            var runner = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)>(), spectra, new List<EngineLayer.DatabaseLoading.DbForTask>(),
                output, new SdrfForRun(Path.Combine(folder, "x.sdrf.tsv"), "provenance"));
            runner.Run();

            Assert.That(File.Exists(Path.Combine(output, SdrfDesignMarker.FileName)), Is.False);
            Assert.That(runner.Warnings, Has.One.StartWith($"Could not copy {SdrfDesignMarker.PathIn(folder)} into the output folder:"));
        }

        /// <summary>
        /// A label-free and an isobaric design both beside the spectra: the run's design check stops on that, and neither is
        /// MetaMorpheus's to choose between, so nothing is written, nothing removed, and nothing refused here.
        /// </summary>
        [Test]
        public static void BothKindsOfDesignBesideTheSpectraAreLeftAlone()
        {
            string folder = NewFolder("BothDesigns");
            var spectra = Spectra(folder, "a1.mzML");
            string sdrf = WriteSdrf(folder, ("a1.raw", "WT", 1));
            string tmtDesign = Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName);
            File.WriteAllText(DesignPath(folder), "label-free");
            File.WriteAllText(tmtDesign, "isobaric");

            var outcome = SdrfAutomaticDesign.Resolve(folder, spectra, SearchTasks(new SearchTask()), new SdrfDesignStatements(), new SdrfDesignWording());

            Assert.That(outcome.Refusal, Is.Null);
            Assert.That(outcome.DesignPath, Is.Null);
            Assert.That(outcome.SdrfPath, Is.EqualTo(Path.GetFullPath(sdrf)));
            Assert.That(File.ReadAllText(DesignPath(folder)), Is.EqualTo("label-free"));
            Assert.That(File.ReadAllText(tmtDesign), Is.EqualTo("isobaric"));
            Assert.That(File.Exists(SdrfDesignMarker.PathIn(folder)), Is.False);
        }

        /// <summary>An SDRF that cannot be read, or has no header, has no comment[plex]; the SDRF reader then says why.</summary>
        [Test]
        public static void AnUnreadableOrEmptySdrfHasNoDefaultPlexColumn()
        {
            string folder = NewFolder("PlexColumnUnread");
            string empty = Path.Combine(folder, "empty.sdrf.tsv");
            File.WriteAllText(empty, string.Empty);

            Assert.That(SdrfAutomaticDesign.DefaultPlexColumnIn(Path.Combine(folder, "missing.sdrf.tsv")), Is.Null);
            Assert.That(SdrfAutomaticDesign.DefaultPlexColumnIn(empty), Is.Null);
            Assert.That(SdrfAutomaticDesign.DefaultPlexColumnIn(WriteTmtSdrf(folder)), Is.EqualTo(SdrfAutomaticDesign.DefaultPlexColumn));
        }

        /// <summary>A marker with any value blank or missing is no marker, so the design beside it is the user's.</summary>
        [TestCase("design")]
        [TestCase("design sha256")]
        [TestCase("sdrf")]
        [TestCase("sdrf sha256")]
        [TestCase("condition")]
        public static void AMarkerMissingAValueIsNoMarker(string key)
        {
            string folder = NewFolder("MarkerMissing_" + key.Replace(' ', '_'));
            new SdrfDesignMarker(GlobalVariables.ExperimentalDesignFileName, "abc", "x.sdrf.tsv", "def", "every factor value column").WriteTo(folder);
            string path = SdrfDesignMarker.PathIn(folder);
            Assume.That(SdrfDesignMarker.ReadFrom(folder), Is.Not.Null);

            File.WriteAllLines(path, File.ReadAllLines(path).Select(l => l.StartsWith(key + "\t") ? key + "\t " : l));
            Assert.That(SdrfDesignMarker.ReadFrom(folder), Is.Null, "blank");

            File.WriteAllLines(path, File.ReadAllLines(path).Where(l => !l.StartsWith(key + "\t")));
            Assert.That(SdrfDesignMarker.ReadFrom(folder), Is.Null, "missing");
        }

        /// <summary>
        /// A marker or design held open by another program cannot be read: the marker is then no marker, and the design is not
        /// known to be unedited, so either way the design is treated as the user's. A design that is gone is not unedited.
        /// </summary>
        [Test]
        public static void AMarkerOrDesignThatCannotBeReadIsNotOurs()
        {
            string folder = NewFolder("Locked");
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1));
            Assert.That(Resolve(folder, spectra, SearchTasks(new SearchTask()), new List<string>()), Is.EqualTo(Continue));
            var marker = SdrfDesignMarker.ReadFrom(folder);
            Assert.That(marker.DescribesUneditedDesignIn(folder), Is.True);

            using (new FileStream(DesignPath(folder), FileMode.Open, FileAccess.Read, FileShare.None))
            {
                Assert.That(marker.DescribesUneditedDesignIn(folder), Is.False, "design locked");
            }
            using (new FileStream(SdrfDesignMarker.PathIn(folder), FileMode.Open, FileAccess.Read, FileShare.None))
            {
                Assert.That(SdrfDesignMarker.ReadFrom(folder), Is.Null, "marker locked");
            }

            File.Delete(DesignPath(folder));
            Assert.That(marker.DescribesUneditedDesignIn(folder), Is.False, "design gone");
        }

        // ---------------------------------------------------------------- through Main

        private static int RunMain(string[] args, out string console)
        {
            TextWriter originalOut = Console.Out;
            var captured = new StringWriter();
            try
            {
                Console.SetOut(captured);
                return Program.Main(args);
            }
            finally
            {
                Console.SetOut(originalOut);
                console = captured.ToString();
            }
        }

        /// <summary>
        /// The wiring: a normalizing run with a refused SDRF beside empty spectra stops with the SDRF's report, before any
        /// search starts, so no output folder is made.
        /// </summary>
        [Test]
        public static void MainStopsANormalizingRunWhoseSdrfGivesNoDesign()
        {
            string folder = NewFolder("MainRefused");
            string spectra = Spectra(folder, "a1.mzML")[0];
            WriteSdrfWithNoFactorColumn(folder, ("a1.raw", 1));
            string database = Path.Combine(folder, "db.fasta");
            File.WriteAllText(database, string.Empty);
            var searchTask = new SearchTask();
            searchTask.SearchParameters.Normalize = true;
            string taskPath = Path.Combine(folder, "SearchTask.toml");
            Toml.WriteFile(searchTask, taskPath, MetaMorpheusTask.tomlConfig);
            string output = Path.Combine(folder, "output");

            int exitCode = RunMain(new[] { "-s", spectra, "-d", database, "-t", taskPath, "-o", output }, out string console);

            Assert.That(exitCode, Is.EqualTo(Stop), console);
            Assert.That(console, Does.Contain("No condition").And.Contain("Normalization needs an experimental design"));
            Assert.That(Directory.Exists(output), Is.False);
        }

        /// <summary>
        /// A whole search with an SDRF beside its spectra: the design and marker are written beside the spectra, the run
        /// quantifies with it, the search writes its own SDRF, results.txt records the input SDRF and how the design was
        /// made, the marker is copied beside Task Settings, and the task file the run saved is the user's, unchanged: the SDRF
        /// reached the search for this run only.
        /// </summary>
        [Test]
        public static void MainSearchesWithTheDesignMadeFromTheSdrf()
        {
            string folder = NewFolder("MainSearch");
            string spectra = Path.Combine(folder, "SmallCalibratible_Yeast.mzML");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "SmallCalibratible_Yeast.mzML"), spectra);
            string sdrf = WriteSdrf(folder, ("SmallCalibratible_Yeast.raw", "WT", 1));
            string database = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "smalldb.fasta");
            string taskPath = Path.Combine(folder, "SearchTask.toml");
            Toml.WriteFile(new SearchTask(), taskPath, MetaMorpheusTask.tomlConfig);
            string output = Path.Combine(folder, "output");

            int exitCode = RunMain(new[] { "-s", spectra, "-d", database, "-t", taskPath, "-o", output, "-v", "minimal" }, out string console);

            Assert.That(exitCode, Is.EqualTo(Continue), console);
            Assert.That(File.Exists(DesignPath(folder)));
            Assert.That(SdrfDesignMarker.ReadFrom(folder).DescribesUneditedDesignIn(folder), Is.True);

            string results = File.ReadAllText(Directory.GetFiles(output, "results.txt", SearchOption.AllDirectories).Single());
            Assert.That(results, Does.Contain($"Input SDRF: {Path.GetFullPath(sdrf)} (sha256 {SdrfDesignMarker.Sha256Of(sdrf)})."));
            Assert.That(results, Does.Contain(SdrfDesignMarker.ReadFrom(folder).Describe()));
            Assert.That(Directory.GetFiles(output, "*.sdrf.tsv", SearchOption.AllDirectories), Has.Length.EqualTo(1),
                "the search's own SDRF, written although the task did not ask for it");

            Assert.That(File.ReadAllText(Path.Combine(output, SdrfDesignMarker.FileName)),
                Is.EqualTo(File.ReadAllText(SdrfDesignMarker.PathIn(folder))), "the marker, beside Task Settings");
            string saved = Directory.GetFiles(Path.Combine(output, "Task Settings"), "*config.toml").Single();
            Assert.That(File.ReadAllText(saved), Does.Not.Contain("RunSdrf").And.Not.Contain(Path.GetFileName(sdrf)));
            var savedTask = MetaMorpheusTask.ReadTaskToml<SearchTask>(saved);
            Assert.That(savedTask.SearchParameters.InputSdrfPath, Is.Null);
            Assert.That(savedTask.SearchParameters.WriteSdrf, Is.False);
        }
    }
}
