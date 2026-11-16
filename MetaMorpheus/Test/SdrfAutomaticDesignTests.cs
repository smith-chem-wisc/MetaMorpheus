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
            Program.ResolveSdrfDesign(folder, spectra, tasks, settings ?? new CommandLineSettings(), normalize, reportToConsole,
                console.Add, () => answer);

        private static string DesignPath(string folder) => Path.Combine(folder, GlobalVariables.ExperimentalDesignFileName);

        // ---------------------------------------------------------------- nothing there, nothing changes

        [Test]
        public static void WithNoSdrfBesideTheSpectraNothingChanges()
        {
            string folder = NewFolder("NoSdrf");
            var spectra = Spectra(folder, "a1.mzML");
            var task = new SearchTask();
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(task), console), Is.EqualTo(Continue));
            Assert.That(Directory.GetFiles(folder).Select(Path.GetFileName), Is.EquivalentTo(new[] { "a1.mzML" }));
            Assert.That(task.SearchParameters.InputSdrfPath, Is.Null);
            Assert.That(task.SearchParameters.WriteSdrf, Is.False);
            Assert.That(console, Is.Empty);
        }

        // ---------------------------------------------------------------- the design and its marker

        /// <summary>
        /// The SDRF names X.raw and the search reads X.mzML: the design names the searched files, MetaMorpheus reads it back
        /// without error, a normalizing run's design check passes, the marker holds both sha256s, the search's own SDRF is
        /// turned on from it, results.txt's line is set, and the SDRF is not touched.
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

            Assert.That(Resolve(folder, spectra, SearchTasks(task), console), Is.EqualTo(Continue), string.Join("\n", console));

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

            Assert.That(task.SearchParameters.InputSdrfPath, Is.EqualTo(Path.GetFullPath(sdrf)));
            Assert.That(task.SearchParameters.WriteSdrf, Is.True);
            Assert.That(task.SearchParameters.DesignWrittenFromSdrf, Is.EqualTo(marker.Describe()));
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
            Assert.That(Resolve(folder, spectra, SearchTasks(task), console), Is.EqualTo(Continue));

            Assert.That(File.ReadAllText(DesignPath(folder)), Is.EqualTo(edited));
            Assert.That(console, Has.Some.Contains("was not written by MetaMorpheus from an SDRF, or was edited since, so it is the design"));
            Assert.That(task.SearchParameters.DesignWrittenFromSdrf, Is.Null);
            Assert.That(task.SearchParameters.InputSdrfPath, Is.Not.Null, "the SDRF still fills the search's own SDRF");
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
                .And.Contain("On the command line, SingleCondition is --sdrfSingleCondition; in a search task's settings, SdrfSingleCondition.")
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

            Assert.That(Resolve(folder, spectra, SearchTasks(task), console, normalize: true), Is.EqualTo(Stop));
            Assert.That(console.Single(), Does.StartWith("2 SDRF files are beside the spectra files (deposit.sdrf.tsv, other.sdrf.tsv)"));
            Assert.That(task.SearchParameters.InputSdrfPath, Is.Null);
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

            Assert.That(Resolve(folder, spectra, SearchTasks(task), console), Is.EqualTo(Continue), string.Join("\n", console));

            Assert.That(File.Exists(DesignPath(folder)), Is.EqualTo(isSdrf));
            Assert.That(task.SearchParameters.InputSdrfPath, Is.EqualTo(isSdrf ? Path.GetFullPath(renamed) : null));
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

            Assert.That(Resolve(folder, spectra, SearchTasks(task), console, normalize: true), Is.EqualTo(Stop));
            Assert.That(console.Single(), Does.StartWith("2 SDRF files are beside the spectra files (deposit.sdrf.tsv, sdrf.tsv)"));
            Assert.That(task.SearchParameters.InputSdrfPath, Is.Null);
            Assert.That(File.Exists(DesignPath(folder)), Is.False);
        }

        // ---------------------------------------------------------------- the grouping the user states

        /// <summary>QP-S42 (b): the task file may state one group, and it is recorded as stated.</summary>
        [Test]
        public static void ATaskSettingStatesOneCondition()
        {
            string folder = NewFolder("TaskSingle");
            var spectra = Spectra(folder, "a1.mzML", "b1.mzML");
            WriteSdrfWithNoFactorColumn(folder, ("a1.raw", 1), ("b1.raw", 2));
            var task = new SearchTask();
            task.SearchParameters.SdrfSingleCondition = "survey";
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(task), console), Is.EqualTo(Continue), string.Join("\n", console));

            var readBack = ExperimentalDesign.ReadExperimentalDesign(DesignPath(folder), spectra, out var errors);
            Assert.That(errors, Is.Empty);
            Assert.That(readBack.Select(f => f.Condition), Is.All.EqualTo("survey"));
            Assert.That(SdrfDesignMarker.ReadFrom(folder).Condition, Is.EqualTo("one condition, 'survey', stated by the user"));
            Assert.That(task.SearchParameters.DesignWrittenFromSdrf, Does.Contain("stated by the user"));
        }

        /// <summary>One design serves every search task, so tasks that group differently are refused, each named.</summary>
        [Test]
        public static void SearchTasksThatGroupDifferentlyAreRefused()
        {
            string folder = NewFolder("TasksDisagree");
            var spectra = Spectra(folder, "a1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1));
            var first = new SearchTask();
            var second = new SearchTask();
            second.SearchParameters.SdrfSingleCondition = "all";
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(first, second), console, normalize: true), Is.EqualTo(Stop));
            Assert.That(console.Last(), Does.StartWith("The search tasks build the SDRF's condition differently (Task1SearchTask: "
                + "every factor value column; Task2SearchTask: one condition, 'all')"));
        }

        /// <summary>The command line's grouping wins over the tasks', and is written into them so their saved settings say so.</summary>
        [Test]
        public static void TheCommandLinesGroupingWinsAndIsRecordedInTheTasks()
        {
            string folder = NewFolder("CommandLineWins");
            var spectra = Spectra(folder, "a1.mzML", "b1.mzML");
            WriteSdrf(folder, ("a1.raw", "WT", 1), ("b1.raw", "KO", 1));
            var task = new SearchTask();
            task.SearchParameters.SdrfSingleCondition = "all";
            var settings = new CommandLineSettings { SdrfConditionColumns = new[] { "genotype" } };
            var console = new List<string>();

            Assert.That(Resolve(folder, spectra, SearchTasks(task), console, settings: settings), Is.EqualTo(Continue), string.Join("\n", console));

            Assert.That(ExperimentalDesign.ReadExperimentalDesign(DesignPath(folder), spectra, out _).Select(f => f.Condition),
                Is.EquivalentTo(new[] { "WT", "KO" }));
            Assert.That(task.SearchParameters.SdrfConditionColumns, Is.EqualTo(new[] { "genotype" }));
            Assert.That(task.SearchParameters.SdrfSingleCondition, Is.Null);
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

            Assert.That(Resolve(folder, spectra, SearchTasks(task), console), Is.EqualTo(Continue));
            Assert.That(task.SearchParameters.InputSdrfPath, Is.EqualTo(own));
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

            Assert.That(Resolve(folder, spectra, SearchTasks(task), console, normalize: true), Is.EqualTo(Continue));
            Assert.That(File.Exists(DesignPath(folder)), Is.False);
            Assert.That(task.SearchParameters.InputSdrfPath, Is.Not.Null);
        }

        // ---------------------------------------------------------------- isobaric

        private static string WriteTmtSdrf(string folder) => WriteSdrf(folder, "deposit.sdrf.tsv",
            new[] { "source name", "characteristics[biological replicate]", "comment[label]", "comment[data file]", "comment[plex]", "factor value[genotype]" },
            new[]
            {
                new[] { "S1", "1", "TMT126", "p1.raw", "P1", "WT" },
                new[] { "S2", "1", "TMT127N", "p1.raw", "P1", "KO" },
            });

        private static SearchTask Tmt11Search()
        {
            var task = new SearchTask();
            task.SearchParameters.DoMultiplexQuantification = true;
            task.SearchParameters.MultiplexModId = "TMT11";
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
        /// quantifies with it, the search writes its own SDRF, and results.txt records the input SDRF and how the design
        /// was made.
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
        }
    }
}
