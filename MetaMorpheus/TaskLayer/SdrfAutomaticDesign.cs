using EngineLayer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using SdrfSearchDesign = Readers.SdrfSearchDesign;
using SdrfSearchDesignOptions = Readers.SdrfSearchDesignOptions;
using MzLibIsobaricMassTag = Omics.Modifications.IsobaricMassTag;

namespace TaskLayer
{
    /// <summary>
    /// What the user stated about an SDRF's design, outside the SDRF: the condition (columns, one condition, or files
    /// whose factor nobody stated as their own condition) and where each file's plex is written. The command line fills
    /// it from its --sdrf... switches. A grouping left empty here is the one stated when MetaMorpheus last wrote the
    /// design, kept in its <see cref="SdrfDesignMarker"/>, or else every factor value column.
    /// </summary>
    public sealed record SdrfDesignStatements
    {
        public IReadOnlyList<string> ConditionColumns { get; init; } = Array.Empty<string>();
        public string SingleCondition { get; init; }
        public bool UnknownFactorIsOwnCondition { get; init; }
        public string PlexColumn { get; init; }
        public string PlexPattern { get; init; }
        public string SinglePlex { get; init; }

        internal bool StatesGrouping => ConditionColumns.Any() || SingleCondition != null || UnknownFactorIsOwnCondition;
    }

    /// <summary>
    /// How a front end names the ways to state a grouping, for the lines <see cref="SdrfAutomaticDesign.Resolve"/> writes:
    /// the command line names its switches, the GUI its own controls.
    /// </summary>
    public sealed record SdrfDesignWording
    {
        /// <summary>Lines added after mzLib's report, saying how to answer a refusal it names (may return empty).</summary>
        public Func<SdrfSearchDesign, string> ReportHint { get; init; } = _ => string.Empty;

        /// <summary>
        /// The end of the line saying a grouping stated for an earlier run is used again: how to state another, or go back
        /// to every factor value column.
        /// </summary>
        public string HowToStateAnotherGrouping { get; init; } =
            "To make it again from every factor value column, delete the design file.";
    }

    /// <summary>
    /// The SDRF beside a run's spectra files, for that run only. The runner (<see cref="EverythingRunnerEngine"/>) hands it
    /// to every search task as <see cref="MetaMorpheusTask.RunSdrf"/>, as it hands them
    /// <see cref="MetaMorpheusTask.AcquiredSpectraFiles"/>: it is never written to a task file, and the task's own settings
    /// stay as the user saved them. A task's own <see cref="SearchParameters.InputSdrfPath"/> still wins.
    /// </summary>
    /// <param name="SdrfPath">The SDRF, as a full path.</param>
    /// <param name="Provenance">
    /// The line results.txt carries when the run's experimental design was written from this SDRF (which SDRF, by sha256,
    /// and how the condition was built); null when it was not.
    /// </param>
    public sealed record SdrfForRun(string SdrfPath, string Provenance)
    {
        /// <summary>The input SDRF a search reads: the task's own when it names one, else the run's SDRF (or null).</summary>
        public static string InputSdrfPathFor(SearchParameters searchParameters, SdrfForRun runSdrf) =>
            !string.IsNullOrWhiteSpace(searchParameters?.InputSdrfPath) ? searchParameters.InputSdrfPath : runSdrf?.SdrfPath;

        /// <summary>
        /// Whether a search writes its own SDRF: when the task asks for it, and whenever the run has an SDRF beside its
        /// spectra files (QuantProject, ruled 2026-11-07).
        /// </summary>
        public static bool WritesSdrfFor(SearchParameters searchParameters, SdrfForRun runSdrf) =>
            searchParameters?.WriteSdrf == true || runSdrf != null;
    }

    /// <summary>What <see cref="SdrfAutomaticDesign.Resolve"/> did, and what a front end has to tell the user.</summary>
    public sealed class SdrfDesignOutcome
    {
        internal SdrfDesignOutcome(List<string> messages, string refusal, string sdrfPath, string designPath,
            string provenance = null)
        {
            Messages = messages;
            Refusal = refusal;
            SdrfPath = sdrfPath;
            DesignPath = designPath;
            RunSdrf = sdrfPath == null ? null : new SdrfForRun(sdrfPath, provenance);
        }

        /// <summary>The lines to show, in order (nothing when no SDRF is beside the spectra files).</summary>
        public IReadOnlyList<string> Messages { get; }

        /// <summary>
        /// Why no design could be made from the SDRF, when one should have been; null otherwise. A front end handles it as
        /// a design file it cannot read: stop when normalizing, otherwise ask whether to continue without a design.
        /// </summary>
        public string Refusal { get; }

        /// <summary>The SDRF beside the spectra files, as a full path; null when there is none (or several).</summary>
        public string SdrfPath { get; }

        /// <summary>The design written from the SDRF; null when none was.</summary>
        public string DesignPath { get; }

        /// <summary>
        /// The SDRF for the run, which a front end passes to <see cref="EverythingRunnerEngine"/>; null when there is none
        /// (or several). Given even when no design was made from it, so the search's own SDRF still copies its samples.
        /// </summary>
        public SdrfForRun RunSdrf { get; }
    }

    /// <summary>
    /// The automatic design (QuantProject gap AUTODESIGN, ruled 2026-11-07, 2026-11-11 and 2026-11-13): when an SDRF is
    /// beside the spectra files (in the first spectra file's folder, named by <see cref="IsSdrfFileName"/>), the run's
    /// experimental design is made from it, and the search's own SDRF is filled from it. With no SDRF there nothing
    /// changes. Shared by the command line and the GUI, which only show the <see cref="SdrfDesignOutcome"/>; it runs before
    /// the run's design check, which then reads the design written here exactly as it reads one the user wrote.
    ///
    /// <list type="bullet">
    /// <item>The SDRF is returned for the run (<see cref="SdrfDesignOutcome.RunSdrf"/>); the front end passes it to the
    /// runner, which hands it to every search task for that run only, so the search's own SDRF is written and carries
    /// every sample fact the input has. The tasks themselves are never changed. A task that names a different SDRF of its
    /// own keeps it, and a line says so.</item>
    /// <item>The design is read with mzLib's SdrfSearchDesign, for exactly the spectra files given. Its condition comes
    /// from <see cref="SdrfDesignStatements"/> when they state one, else from the grouping stated when MetaMorpheus last
    /// wrote this design (kept in its marker, and said when used), else from every factor value column. The kit is the
    /// search tasks' multiplex label; the plex is read from a stated column, pattern or single plex, else from the SDRF's
    /// <see cref="DefaultPlexColumn"/>.</item>
    /// <item>The design is written beside the first spectra file with a <see cref="SdrfDesignMarker"/> holding its and
    /// the SDRF's sha256 (QP-S41). A design already there that the marker shows MetaMorpheus wrote, unedited, is written
    /// again for this run's files; any other design there is the user's, is used as is and never touched, and a line says
    /// the SDRF was not read for it.</item>
    /// <item>A design that cannot be made (refused, unreadable, two SDRFs, search tasks of different kits, or a design or
    /// marker that cannot be written) is returned as a <see cref="SdrfDesignOutcome.Refusal"/>. Nothing is written then,
    /// and a design MetaMorpheus wrote before, which no longer fits, is removed with its marker; a grouping the marker kept
    /// goes with it, and a line says so. The SDRF is never edited.</item>
    /// <item>Two SDRFs are never chosen between, but they are refused only where a design would be made from one: a search
    /// that quantifies nothing, and a design that is the user's, go on as with no SDRF.</item>
    /// </list>
    /// </summary>
    public static class SdrfAutomaticDesign
    {
        /// <summary>The ending every SDRF file name carries (the SDRF-Proteomics specification's).</summary>
        public const string SdrfFileEnding = ".sdrf.tsv";

        /// <summary>
        /// A file named just this is an SDRF too (QuantProject, ruled 2026-11-13). It lacks the dot
        /// <see cref="SdrfFileEnding"/> starts with, and it is the name a person most often gives one, so missing it would
        /// make no design and say nothing.
        /// </summary>
        public const string SdrfBareFileName = "sdrf.tsv";

        /// <summary>
        /// TMT16 has no multiplex label of its own in MetaMorpheus: its channels are TMT18's first sixteen, so the file is
        /// one a TMT18 run reads as is. Said, so the search is not set to a label that is not there.
        /// </summary>
        public const string Tmt16Note = "TMT16 is searched in MetaMorpheus with its TMT18 multiplex label (TMTpro, the same reagent). "
            + "This design is read as is by such a run; channels 134C and 135N stay empty.";

        /// <summary>
        /// The column a TMT SDRF's plex is read from when no plex source is stated (QuantProject, ruled 2026-11-13): the
        /// column PXReprise writes in every TMT SDRF it drafts. A fixed column name is the SDRF stating the plex, not
        /// MetaMorpheus guessing it; a stated column, file-name pattern or single plex wins over it.
        /// </summary>
        public const string DefaultPlexColumn = "comment[plex]";

        /// <summary>
        /// The SDRF's own spelling of <see cref="DefaultPlexColumn"/> (case ignored), or null when it has no such column or
        /// cannot be read (the SDRF reader then says why).
        /// </summary>
        public static string DefaultPlexColumnIn(string sdrfPath)
        {
            string header;
            try
            {
                header = File.ReadLines(sdrfPath).FirstOrDefault() ?? string.Empty;
            }
            catch (Exception e) when (e is IOException or UnauthorizedAccessException or ArgumentException or NotSupportedException)
            {
                return null;
            }
            return header.TrimStart('\uFEFF').Split('\t').Select(c => c.Trim())
                .FirstOrDefault(c => string.Equals(c, DefaultPlexColumn, StringComparison.OrdinalIgnoreCase));
        }

        /// <summary>
        /// Where a TMT design's plex is read: the stated column, pattern or single plex when one is stated, otherwise the
        /// SDRF's <see cref="DefaultPlexColumn"/> when it has one, otherwise nothing (mzLib then refuses, naming what is
        /// missing; <see cref="DefaultPlexNote"/> adds that comment[plex] was looked for).
        /// </summary>
        public static (string Column, string Pattern, string SinglePlex) PlexSource(string sdrfPath, string column, string pattern,
            string singlePlex) =>
            PlexIsStated(column, pattern, singlePlex)
                ? (column, pattern, singlePlex)
                : (DefaultPlexColumnIn(sdrfPath), null, null);

        /// <summary>True when the user stated where the plex is read: a column, a file-name pattern or a single plex.</summary>
        public static bool PlexIsStated(string column, string pattern, string singlePlex) =>
            column != null || pattern != null || singlePlex != null;

        /// <summary>
        /// A line saying where an isobaric design's plex came from when the user stated none (empty otherwise). mzLib's report
        /// calls any plex column it is given "declared", so without this line a plex MetaMorpheus chose reads as the user's.
        /// When the SDRF has no comment[plex], mzLib refuses (no plex source was declared), and the line says comment[plex]
        /// was looked for; otherwise it says comment[plex] was read.
        /// </summary>
        public static string DefaultPlexNote(SdrfSearchDesign design, bool isobaric, bool plexStated)
        {
            if (!isobaric || plexStated)
                return string.Empty;
            return design.Refusals.Any(r => r.StartsWith("No plex source was declared", StringComparison.Ordinal))
                ? $"The SDRF has no {DefaultPlexColumn} column, which MetaMorpheus reads for the plex when none is named."
                    + Environment.NewLine
                : $"No plex source was named, so MetaMorpheus read the plex from the SDRF's {DefaultPlexColumn} column."
                    + Environment.NewLine;
        }

        /// <summary>
        /// True for a file the automatic design reads as an SDRF: a name ending <see cref="SdrfFileEnding"/>, or exactly
        /// <see cref="SdrfBareFileName"/>. Case is ignored.
        /// </summary>
        public static bool IsSdrfFileName(string path)
        {
            string name = Path.GetFileName(path);
            // macOS leaves an AppleDouble file named ._<name> beside a file it copies to some drives and shares; it is
            // metadata, not an SDRF (pcruzparri's review of #2980).
            if (name.StartsWith("._", StringComparison.Ordinal))
                return false;
            return name.EndsWith(SdrfFileEnding, StringComparison.OrdinalIgnoreCase)
                || string.Equals(name, SdrfBareFileName, StringComparison.OrdinalIgnoreCase);
        }

        /// <summary>
        /// Makes the run's design from the SDRF beside the spectra files, when there is one (see the class summary), and says
        /// what it did. Changes only the design files beside the spectra; never the search tasks or the SDRF.
        /// </summary>
        /// <param name="designDirectory">The folder of the first spectra file, where a run reads its design.</param>
        /// <param name="spectraFilePaths">Full paths of the spectra files this run will search.</param>
        /// <param name="searchTasks">The run's search tasks, by name; the design serves them all.</param>
        /// <param name="stated">What the user stated outside the SDRF and the tasks (the command line's switches).</param>
        /// <param name="wording">How the front end names the ways to state a grouping.</param>
        public static SdrfDesignOutcome Resolve(string designDirectory, List<string> spectraFilePaths,
            IReadOnlyList<(string Name, SearchTask Task)> searchTasks, SdrfDesignStatements stated, SdrfDesignWording wording)
        {
            stated ??= new SdrfDesignStatements();
            wording ??= new SdrfDesignWording();
            var messages = new List<string>();

            if (searchTasks.Count == 0)
            {
                return new SdrfDesignOutcome(messages, null, null, null);
            }

            var sdrfs = Directory.GetFiles(designDirectory)
                .Where(IsSdrfFileName)
                .OrderBy(f => f, StringComparer.OrdinalIgnoreCase)
                .ToList();
            if (sdrfs.Count == 0)
            {
                return new SdrfDesignOutcome(messages, null, null, null);
            }

            // Which SDRF, when a folder holds several, is not ruled (QuantProject, 2026-11-07): one is kept per group of
            // files. So none is guessed, and none goes to the search; they are refused below only where a design would be
            // made from one (pcruzparri's review of #2980).
            string sdrf = sdrfs.Count == 1 ? Path.GetFullPath(sdrfs[0]) : null;
            string several = sdrf == null
                ? $"{sdrfs.Count} SDRF files are beside the spectra files ({string.Join(", ", sdrfs.Select(Path.GetFileName))}), "
                    + "and MetaMorpheus does not choose between them"
                : null;
            if (sdrf != null)
            {
                NoteSdrfForSearchTasks(sdrf, searchTasks, messages);
            }

            // A search that quantifies nothing reads no design; the kit is read from the tasks that quantify.
            var quantifying = searchTasks
                .Where(t => t.Task.SearchParameters.DoLabelFreeQuantification || t.Task.SearchParameters.DoMultiplexQuantification)
                .ToList();
            if (quantifying.Count == 0)
            {
                if (several != null)
                {
                    messages.Add(several + ", so neither is given to the search.");
                }
                return new SdrfDesignOutcome(messages, null, sdrf, null);
            }

            string labelFreePath = Path.Combine(designDirectory, GlobalVariables.ExperimentalDesignFileName);
            string isobaricPath = Path.Combine(designDirectory, GlobalVariables.TmtExperimentalDesignFileName);
            string existing = new[] { labelFreePath, isobaricPath }.Where(File.Exists).ToList() switch
            {
                { Count: 0 } => null,
                { Count: 1 } one => one[0],
                // both kinds: the design check stops the run, and neither is ours to choose between
                _ => string.Empty,
            };
            if (existing == string.Empty)
            {
                return new SdrfDesignOutcome(messages, null, sdrf, null);
            }

            var marker = SdrfDesignMarker.ReadFrom(designDirectory);
            bool ours = existing != null && marker != null
                && string.Equals(marker.DesignFileName, Path.GetFileName(existing), StringComparison.Ordinal)
                && marker.DescribesUneditedDesignIn(designDirectory);
            if (existing != null && !ours)
            {
                messages.Add($"{Path.GetFileName(existing)} is beside the spectra files and was not written by MetaMorpheus from an SDRF, "
                    + "or was edited since, so it is the design; "
                    + (sdrf != null ? $"the SDRF beside it ({sdrf}) is not read for the design." : $"the {sdrfs.Count} SDRF files beside it are not read."));
                return new SdrfDesignOutcome(messages, null, sdrf, null);
            }

            string removeOurs = ours ? existing : null;
            var ourMarker = ours ? marker : null;

            if (several != null)
            {
                return Refuse(messages, several + ", so no experimental design was made from an SDRF. Keep beside the spectra files "
                    + "only the SDRF that describes them.", null, null, removeOurs, designDirectory, ourMarker);
            }

            var grouping = ReadGrouping(stated, ourMarker, wording, messages);

            var kits = quantifying
                .Select(t => t.Task.SearchParameters.DoMultiplexQuantification ? t.Task.SearchParameters.MultiplexModId ?? string.Empty : null)
                .Distinct()
                .ToList();
            if (kits.Count > 1)
            {
                return Refuse(messages,
                    "The search tasks quantify differently (" + string.Join(", ", quantifying.Select(t => t.Name + ": "
                        + (t.Task.SearchParameters.DoMultiplexQuantification ? t.Task.SearchParameters.MultiplexModId : "label-free")))
                    + "), and one experimental design cannot serve them all, so none was made from the SDRF.",
                    null, sdrf, removeOurs, designDirectory, ourMarker);
            }

            MzLibIsobaricMassTag tag = null;
            if (kits[0] != null && !MzLibIsobaricMassTag.TryGetIsobaricMassTag(kits[0], GlobalVariables.AllModsKnown, out tag))
            {
                return Refuse(messages,
                    $"The search's multiplex label '{kits[0]}' is not an isobaric kit mzLib knows, so the SDRF's channels cannot be read.",
                    null, sdrf, removeOurs, designDirectory, ourMarker);
            }

            var plex = tag != null
                ? PlexSource(sdrf, stated.PlexColumn, stated.PlexPattern, stated.SinglePlex)
                : (stated.PlexColumn, stated.PlexPattern, stated.SinglePlex);
            var options = new SdrfSearchDesignOptions
            {
                Tag = tag,
                PlexColumn = plex.Item1,
                PlexFileNamePattern = plex.Item2,
                SinglePlex = plex.Item3,
                ConditionColumns = grouping.Columns,
                SingleCondition = grouping.SingleCondition,
                UnknownFactorIsOwnCondition = grouping.UnknownIsOwnCondition,
            };
            var design = SdrfSearchDesign.Read(sdrf, spectraFilePaths, options);
            string report = design.Report()
                + DefaultPlexNote(design, tag != null, PlexIsStated(stated.PlexColumn, stated.PlexPattern, stated.SinglePlex))
                + wording.ReportHint(design);

            if (!design.IsValid)
            {
                return Refuse(messages,
                    $"An SDRF is beside the spectra files ({sdrf}), but no experimental design could be made from it, for the reasons above.",
                    report, sdrf, removeOurs, designDirectory, ourMarker);
            }

            string designPath = tag != null ? isobaricPath : labelFreePath;
            SdrfDesignMarker newMarker;
            try
            {
                if (existing != null && existing != designPath)
                {
                    File.Delete(existing); // ours, unedited, and of the other kind: the search's label changed
                }
                design.WriteDesign(designPath);

                newMarker = new SdrfDesignMarker(Path.GetFileName(designPath), SdrfDesignMarker.ContentSha256Of(designPath), sdrf,
                    SdrfDesignMarker.Sha256Of(sdrf), DescribeCondition(design, grouping),
                    grouping.Columns, grouping.SingleCondition, grouping.UnknownIsOwnCondition);
                newMarker.WriteTo(designDirectory);
            }
            catch (Exception e) when (e is IOException or UnauthorizedAccessException)
            {
                // A spectra folder or design MetaMorpheus cannot write: refused like any design it cannot make, never a crash
                // (pcruzparri's review of #2980). Each file here is MetaMorpheus's: the user's design returned above.
                messages.Add(report.TrimEnd());
                messages.AddRange(RemoveOurs(new[] { designPath, existing }, designDirectory));
                return new SdrfDesignOutcome(messages,
                    $"An SDRF is beside the spectra files ({sdrf}), but the experimental design made from it could not be written there: "
                    + e.Message, sdrf, null);
            }

            messages.Add(report.TrimEnd());
            if (tag?.TagType == Omics.Modifications.IsobaricMassTagType.TMT16)
            {
                messages.Add(Tmt16Note);
            }
            messages.Add($"Wrote {designPath} from the SDRF {sdrf}" + (ours ? ", again for this run's spectra files." : "."));
            return new SdrfDesignOutcome(messages, null, sdrf, designPath, newMarker.Describe());
        }

        /// <summary>
        /// A design that could not be made from the SDRF. A design MetaMorpheus wrote from the SDRF before, unedited, is
        /// removed with its marker first: it no longer fits, and the run must not quantify with it.
        /// </summary>
        private static SdrfDesignOutcome Refuse(List<string> messages, string why, string report, string sdrf, string removeOurs,
            string designDirectory, SdrfDesignMarker ourMarker)
        {
            if (report != null)
            {
                messages.Add(report.TrimEnd());
            }
            if (removeOurs != null)
            {
                var notRemoved = RemoveOurs(new[] { removeOurs }, designDirectory);
                if (notRemoved.Count == 0)
                {
                    // The marker was the only place a stated grouping was kept, so the next run will not use it again unless
                    // it is stated again (pcruzparri's review of #2980).
                    messages.Add($"{removeOurs} was removed: MetaMorpheus wrote it from the SDRF before, and it no longer fits."
                        + (ourMarker?.StatesGrouping == true
                            ? $" The grouping kept with it ({DescribeGrouping(new SdrfGrouping(ourMarker.ConditionColumns.ToList(), ourMarker.SingleCondition, ourMarker.UnknownFactorIsOwnCondition))}) "
                                + "went with its marker; state it again for the next run."
                            : string.Empty));
                }
                else
                {
                    messages.AddRange(notRemoved);
                }
            }
            return new SdrfDesignOutcome(messages, why, sdrf, null);
        }

        /// <summary>
        /// Removes designs MetaMorpheus wrote beside the spectra, then their marker, and returns a line for each file that
        /// could not be removed (empty when all are gone). The marker goes only once every design has: a design left
        /// without its marker would be taken for the user's by the next run and used as is, while with it the next run
        /// knows the design is MetaMorpheus's and makes it again or refuses again.
        /// </summary>
        private static List<string> RemoveOurs(IEnumerable<string> designs, string designDirectory)
        {
            var notRemoved = designs
                .Where(p => p != null)
                .Distinct(StringComparer.OrdinalIgnoreCase)
                .Select(TryDelete)
                .Where(problem => problem != null)
                .Select(problem => problem + " A run that goes on reads it as the design, so remove it by hand.")
                .ToList();
            if (notRemoved.Count == 0)
            {
                string markerProblem = TryDelete(SdrfDesignMarker.PathIn(designDirectory));
                if (markerProblem != null)
                {
                    notRemoved.Add(markerProblem);
                }
            }
            return notRemoved;
        }

        /// <summary>
        /// Deletes a file MetaMorpheus wrote beside the spectra (a file already gone is fine): null when it is gone, otherwise
        /// the line saying it could not be removed, and why.
        /// </summary>
        private static string TryDelete(string path)
        {
            try
            {
                File.Delete(path);
                return null;
            }
            catch (Exception e) when (e is IOException or UnauthorizedAccessException)
            {
                return $"{path} could not be removed: {e.Message}";
            }
        }

        /// <summary>The grouping the user states for a design made from an SDRF: the columns, one condition, or neither.</summary>
        private sealed record SdrfGrouping(List<string> Columns, string SingleCondition, bool UnknownIsOwnCondition);

        /// <summary>
        /// The grouping for the design: the one stated now (the command line's switches) when there is one; otherwise the one
        /// stated when MetaMorpheus last wrote this design, kept in its <paramref name="ourMarker"/> (a line says it is used
        /// again, and how to state another); otherwise every factor value column. One design serves the whole run, so the
        /// grouping belongs with the design, never with a search task (nbollis's review of #2980).
        /// </summary>
        /// <param name="ourMarker">The marker of a design MetaMorpheus wrote and nobody edited since; null otherwise.</param>
        private static SdrfGrouping ReadGrouping(SdrfDesignStatements stated, SdrfDesignMarker ourMarker, SdrfDesignWording wording,
            List<string> messages)
        {
            if (stated.StatesGrouping)
            {
                return new SdrfGrouping(stated.ConditionColumns.ToList(), stated.SingleCondition, stated.UnknownFactorIsOwnCondition);
            }

            if (ourMarker?.StatesGrouping == true)
            {
                var kept = new SdrfGrouping(ourMarker.ConditionColumns.ToList(), ourMarker.SingleCondition, ourMarker.UnknownFactorIsOwnCondition);
                messages.Add($"The grouping stated when MetaMorpheus last wrote {ourMarker.DesignFileName} from the SDRF is used again "
                    + $"({DescribeGrouping(kept)}), as {SdrfDesignMarker.FileName} records it. {wording.HowToStateAnotherGrouping}");
                return kept;
            }

            return new SdrfGrouping(new List<string>(), null, false);
        }

        private static string DescribeGrouping(SdrfGrouping g) =>
            (g.SingleCondition != null ? $"one condition, '{g.SingleCondition}'"
                : g.Columns.Count > 0 ? "columns " + string.Join(", ", g.Columns) : "every factor value column")
            + (g.UnknownIsOwnCondition ? ", unknown factors their own condition" : string.Empty);

        /// <summary>How the condition of a design was built, for its marker and results.txt: read or stated.</summary>
        private static string DescribeCondition(SdrfSearchDesign design, SdrfGrouping grouping)
        {
            string single = design.LabelFree?.SingleCondition ?? design.Isobaric?.SingleCondition;
            var columns = design.LabelFree?.ConditionColumns ?? design.Isobaric?.ConditionColumns ?? Array.Empty<string>();
            string text = single != null
                ? $"one condition, '{single}', stated by the user"
                : (grouping.Columns.Count > 0 ? "the column(s) the user named, " : "every factor value column, ") + string.Join(" + ", columns);
            return grouping.UnknownIsOwnCondition
                ? text + "; files whose factor nobody stated are a condition of their own, as stated by the user"
                : text;
        }

        /// <summary>
        /// Says what the SDRF beside the spectra files does to each search task in this run, which the runner hands it to
        /// (<see cref="SdrfForRun"/>): the search's own SDRF is written whether or not the task asked for it (QuantProject,
        /// ruled 2026-11-07), and a task naming another SDRF of its own keeps it. The tasks are not changed.
        /// </summary>
        private static void NoteSdrfForSearchTasks(string sdrf, IReadOnlyList<(string Name, SearchTask Task)> searchTasks,
            List<string> messages)
        {
            foreach (var (name, task) in searchTasks)
            {
                string own = task.SearchParameters.InputSdrfPath;
                if (!string.IsNullOrWhiteSpace(own)
                    && !string.Equals(Path.GetFullPath(own.Trim()), sdrf, StringComparison.OrdinalIgnoreCase))
                {
                    messages.Add($"{name} names its own InputSdrfPath ({Path.GetFullPath(own.Trim())}), not the SDRF beside the spectra "
                        + $"files ({sdrf}); the search's SDRF copies its samples from the task's own.");
                }

                if (!task.SearchParameters.WriteSdrf)
                {
                    messages.Add($"{name} writes the search's own SDRF, because an SDRF is beside the spectra files.");
                }
            }
        }
    }
}
