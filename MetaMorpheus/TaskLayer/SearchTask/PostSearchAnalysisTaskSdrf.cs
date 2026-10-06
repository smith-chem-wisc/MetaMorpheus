using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using Omics.Modifications;
using Proteomics;
using Readers;
using UsefulProteomicsDatabases;
using IsobaricMassTag = EngineLayer.IsobaricMassTag;
using IsobaricMassTagType = EngineLayer.IsobaricMassTagType;

namespace TaskLayer
{
    /// <summary>
    /// Writes an SDRF-Proteomics file describing the search, next to its results.
    ///
    /// All of the format, vocabulary and validation logic lives in mzLib
    /// (<see cref="SdrfBuilder"/>); this is the adapter that hands it what MetaMorpheus knows.
    /// Deliberately thin: anything that reasons about SDRF belongs upstream, where it can be tested
    /// against the curated corpus rather than through a search.
    /// </summary>
    public partial class PostSearchAnalysisTask
    {
        /// <summary>
        /// Written per search, into that search's own output folder, and never merged in place.
        /// Repeated searches of one experiment therefore accumulate as separate files that pool by
        /// accession rather than overwriting each other -- and two concurrent searches over the same
        /// spectra cannot race on a shared file.
        /// </summary>
        private void WriteSdrf()
        {
            // Deliberately wrapped, unlike the other writers in this class, which are bare and abort
            // the whole run on failure. An SDRF is metadata about results, not results: a bad sample
            // annotation must not destroy a search that has already succeeded. Missing metadata was
            // already named before the run, by SearchTask.WarnAboutSdrfGaps.
            try
            {
                var rows = BuildSdrfRows().ToList();
                if (rows.Count == 0)
                {
                    Warn("No spectra files to describe; skipping SDRF.");
                    return;
                }

                var document = SdrfBuilder.Build(rows, new SdrfBuilderOptions
                {
                    Software = new CvParam("MS", "MS:1002826", "MetaMorpheus", ""),
                    SoftwareVersion = GlobalVariables.MetaMorpheusVersion,
                    // The sample metadata is only as complete as the input allowed. Gaps were named
                    // before the run; by here it is committed, so accept what we have and let
                    // SdrfCoverage report on it rather than throwing mid-write.
                    RequireSampleMetadata = false,
                    // comment[label] is written bare (TMT126), not accessioned: it is what the
                    // community writes, and quantms crashes on the accessioned form.
                    LabelForm = SdrfLabelForm.Bare
                });

                string path = Path.Combine(Parameters.OutputFolder, "experiment.sdrf.tsv");
                document.WriteResults(path);
                FinishedWritingFile(path, new List<string> { Parameters.SearchTaskId });

                ReportSdrfCoverage(document);
            }
            catch (Exception e)
            {
                EngineCrashed("SdrfWriter", e);
            }
        }

        /// <summary>
        /// One row per spectra file, pairing the sample facts with the assay facts.
        /// </summary>
        private IEnumerable<SdrfRowInput> BuildSdrfRows()
        {
            // The sample half. ExperimentalDesign.tsv is the only place MetaMorpheus holds it, and
            // it is OPTIONAL: when absent, PostSearchAnalysisTask fabricates a degenerate design in
            // which every file is its own biological replicate with an empty condition. That design is
            // not read here. Without a design the replicate and fraction numbers are UNKNOWN, so they
            // are passed as null and written `not available` (mzLib #1378), never a claimed 1.
            //
            // An isobaric search keeps its design in TmtDesign.txt instead, and that file is the only
            // place the channel-to-sample map exists. When it is usable the SDRF gets one row per
            // sample per channel, as the specification wants; when it is not, the search is described
            // one row per file with the label left unresolved, exactly as before. A file whose channel
            // rows could not say which channel is which (ChannelRowsUnusableReason) gets the same
            // fallback, keeping the fraction and technical replicate the design gives it.
            //
            // ExperimentalDesign.tsv is never read for an isobaric search, even when TmtDesign.txt is
            // missing or unusable: SearchTask.WarnAboutSdrfGaps told the user before the run that it is
            // not consulted, and a file-level condition would misdescribe a file holding many samples.
            var isobaric = ReadIsobaricDesignIfPresent(out var tagType);
            var design = Parameters.SearchParameters.DoMultiplexQuantification
                ? new Dictionary<string, SpectraFileInfo>(StringComparer.OrdinalIgnoreCase)
                : ReadExperimentalDesignIfPresent();
            var describedOnce = new List<string>();

            // A sample name reused across plexes for different material (#2817 review): only those names
            // are scoped to their plex. A bridge, which agrees on condition and replicate, keeps its name.
            var reusedSampleNames = isobaric is null
                ? new HashSet<string>()
                : SampleNamesReusedForDifferentSamples(isobaric.Values
                    .Where(f => ChannelRowsUnusableReason(f, tagType!.Value) is null));
            if (reusedSampleNames.Any())
                Warn(ReusedSampleNamesWarning(reusedSampleNames, SampleNamesRepeatedWithinAPlex(isobaric!.Values
                    .Where(f => ChannelRowsUnusableReason(f, tagType!.Value) is null))));

            var organism = ResolveOrganismFromSearchDatabase();

            foreach (var rawFilePath in Parameters.CurrentRawFileList)
            {
                // comment[data file] is the ACQUIRED file and comment[searched data file] the derivative a
                // Calibrate or Average task made (sdrf D46); BuildAssay writes both, for file and channel rows.
                string stem = Path.GetFileNameWithoutExtension(rawFilePath);

                // Per-file parameters, not task-level. MetaMorpheus supports per-file overrides of
                // protease, tolerances and dissociation type, and using the task-level values would
                // silently flatten real differences between rows.
                //
                // FileName here is the FULL PATH the task was given (MetaMorpheusTask stores
                // currentRawDataFilepathList[i]), not a bare file name, and it comes from the same
                // list being walked -- so match on rawFilePath. Matching on the bare name never
                // hit, and every row quietly reported the task-level values instead.
                var common = FileSpecificParameters
                    ?.FirstOrDefault(f => string.Equals(f.FileName, rawFilePath, StringComparison.OrdinalIgnoreCase))
                    .Parameters ?? CommonParameters;

                TmtFileInfo tmtFile = null;
                if (isobaric is not null && isobaric.TryGetValue(Path.GetFullPath(rawFilePath), out tmtFile))
                {
                    if (ChannelRowsUnusableReason(tmtFile, tagType!.Value) is { } reason)
                    {
                        describedOnce.Add(reason);
                    }
                    else
                    {
                        foreach (var row in BuildChannelRows(tmtFile, tagType.Value, organism,
                                     BuildAssay(rawFilePath, common, tmtFile.TechnicalReplicate, tmtFile.Fraction),
                                     reusedSampleNames))
                            yield return row;
                        continue;
                    }
                }

                design.TryGetValue(stem, out var sampleInfo);

                var sample = new SdrfSample
                {
                    // Written, as `not available` when unknown: the specification requires both columns, and a
                    // column that is present says "nobody filled this in" where an absent one says nothing.
                    Characteristics = new Dictionary<string, CvParam>
                    {
                        ["characteristics[disease]"] = null,
                        ["characteristics[cell type]"] = null
                    },
                    SourceName = sampleInfo?.Condition is { Length: > 0 } condition
                        ? condition + " " + (sampleInfo.BiologicalReplicate + 1)
                        : stem,
                    Organism = organism,
                    // SDRF is 1-based; SpectraFileInfo stores these 0-based. With no design nobody established
                    // the numbers, so null, which mzLib writes as `not available`, rather than a claimed 1.
                    BiologicalReplicate = sampleInfo?.BiologicalReplicate + 1,
                    Label = ResolveLabel(Parameters.SearchParameters),
                    FactorValue = sampleInfo?.Condition,
                    FactorValueColumn = string.IsNullOrWhiteSpace(sampleInfo?.Condition)
                        ? null
                        : "factor value[condition]"
                };

                var assay = BuildAssay(rawFilePath, common,
                    technicalReplicate: tmtFile?.TechnicalReplicate ?? sampleInfo?.TechnicalReplicate + 1,
                    fraction: tmtFile?.Fraction ?? sampleInfo?.Fraction + 1);

                yield return new SdrfRowInput(sample, assay);
            }

            if (describedOnce.Any())
                Warn("The SDRF describes these files once, without their channels or samples: " +
                     string.Join("; ", describedOnce.Distinct()) + ".");
        }

        /// <summary>
        /// Why channel rows cannot be written for <paramref name="file"/>, or null when they can. Each
        /// case falls back to one row per file, because a channel row that cannot say which channel it
        /// is misdescribes the file:
        ///
        /// 1. A tag type PRIDE has no channel terms for (DiLeu). Every row would carry
        ///    `not available` in comment[label], which mzLib's SdrfQuantAuditor reads as a label-free
        ///    file holding that many samples.
        /// 2. A plex with no annotated channels: the placeholder row TmtExperimentalDesign.Write emits,
        ///    which Read accepts without error. There are no channel rows to write, and the file used to
        ///    vanish from the SDRF.
        /// 3. A channel the search's plex does not have (a TMT6 "127" on a TMT11 search). Read does not
        ///    know the tag type; this is the check <see cref="TmtExperimentalDesign.ToMzLibDesign"/>
        ///    makes, and quantification skips such a channel. Written, it would carry another kit's label.
        ///
        /// Shared with SearchTask's pre-run warning, so a search is told before it runs.
        /// </summary>
        internal static string ChannelRowsUnusableReason(TmtFileInfo file, IsobaricMassTagType tagType)
        {
            if (PrideLabelFamily(tagType) is null)
                return $"the PRIDE vocabulary defines channel terms for TMT and iTRAQ only, not {tagType}, " +
                       "so comment[label] could not say which channel a row is";

            string fileName = Path.GetFileName(file.FullFilePathWithExtension);
            if (file.Annotations is not { Count: > 0 })
                return $"{fileName} has no annotated channels in {GlobalVariables.TmtExperimentalDesignFileName}";

            var channels = IsobaricMassTag.GetReporterIonLabels(tagType) ?? new List<string>();
            var offPlex = file.Annotations
                .Select(a => a.Tag?.Trim() ?? "")
                .Where(tag => !channels.Contains(tag, StringComparer.OrdinalIgnoreCase))
                .ToList();
            return offPlex.Any()
                ? $"{fileName} annotates {string.Join(", ", offPlex.Select(t => "'" + t + "'"))}, which is not a " +
                  $"channel of {tagType}"
                : null;
        }

        /// <summary>
        /// The channel sample names that name more than one sample: a name written in more than one plex,
        /// or on more than one channel of one plex, with a different (condition, biological replicate) on each. <see cref="TmtExperimentalDesign.Read"/>
        /// checks names only within a plex, on purpose -- a bridge channel carries one name across plexes --
        /// so a design that names each plex's channels S1, S2, ... for different material loads cleanly.
        /// Written as is, one source name would claim two conditions, and anything keyed on source name
        /// (REQ-2, MAP-01) would merge two samples. A name that agrees everywhere is a bridge and is not
        /// returned. Blank names are not returned: they are already scoped to their plex.
        ///
        /// Shared with SearchTask's pre-run warning, so a search is told before it runs.
        /// </summary>
        internal static HashSet<string> SampleNamesReusedForDifferentSamples(IEnumerable<TmtFileInfo> files) =>
            files
                .SelectMany(f => f.Annotations)
                .Where(a => !string.IsNullOrWhiteSpace(a.SampleName))
                .GroupBy(a => a.SampleName, StringComparer.Ordinal)
                .Where(g => g.Select(a => ((a.Condition ?? "").Trim(), a.BiologicalReplicate)).Distinct().Count() > 1)
                .Select(g => g.Key)
                .ToHashSet(StringComparer.Ordinal);

        /// <summary>
        /// The sample names one plex gives to more than one channel, with a different (condition, biological
        /// replicate) on each: what Auto-fill writes when every channel is named after the cell line.
        /// <see cref="TmtExperimentalDesign.Read"/> checks uniqueness within a plex on (sample, biological
        /// replicate, fraction, technical replicate), not on the name, so this loads cleanly. Adding the plex
        /// alone would leave these channels one source name, so they are scoped to the channel as well
        /// (pcruzparri, #2817 review). Every fraction of a plex carries the same annotations, so the scope
        /// is stable across fractions.
        /// </summary>
        internal static HashSet<string> SampleNamesRepeatedWithinAPlex(IEnumerable<TmtFileInfo> files) =>
            files
                .SelectMany(f => RepeatedWithin(f.Annotations))
                .ToHashSet(StringComparer.Ordinal);

        private static IEnumerable<string> RepeatedWithin(IEnumerable<TmtPlexAnnotation> annotations) =>
            annotations
                .Where(a => !string.IsNullOrWhiteSpace(a.SampleName))
                .GroupBy(a => a.SampleName, StringComparer.Ordinal)
                .Where(g => g.Select(a => ((a.Condition ?? "").Trim(), a.BiologicalReplicate)).Distinct().Count() > 1)
                .Select(g => g.Key);

        /// <summary>
        /// The warning for names that name more than one sample. A name repeated inside one plex is
        /// described as that, not as reuse across plexes, because the fix the user needs is different.
        /// </summary>
        internal static string ReusedSampleNamesWarning(IEnumerable<string> names, ISet<string> repeatedWithinAPlex)
        {
            string Quote(IEnumerable<string> list) =>
                string.Join(", ", list.OrderBy(n => n, StringComparer.Ordinal).Select(n => "'" + n + "'"));

            var withinPlex = names.Where(n => repeatedWithinAPlex?.Contains(n) == true).ToList();
            var acrossPlexes = names.Except(withinPlex, StringComparer.Ordinal).ToList();
            var parts = new List<string>();
            if (acrossPlexes.Any())
                parts.Add("The design uses these sample names in more than one plex with a different condition or " +
                          "biological replicate, so they name different samples: " + Quote(acrossPlexes) +
                          ". The SDRF writes each as '<plex> <sample name>' so that one source name stays one sample.");
            if (withinPlex.Any())
                parts.Add("The design gives these sample names to more than one channel of the same plex, with a " +
                          "different condition or biological replicate on each, so they name different samples: " +
                          Quote(withinPlex) + ". The SDRF writes each as '<plex> <sample name> <channel>' so that " +
                          "one source name stays one sample.");
            parts.Add("Give them distinct names in " + GlobalVariables.TmtExperimentalDesignFileName +
                      " to choose the names yourself.");
            return string.Join(" ", parts);
        }

        /// <summary>
        /// The assay half of a row: everything the search itself knows about one data file. Shared by
        /// every channel row of an isobaric file, which is what makes those rows one assay.
        /// </summary>
        private SdrfAssay BuildAssay(string rawFilePath, CommonParameters common, int? technicalReplicate, int? fraction) =>
            new SdrfAssay
            {
                // Two names (sdrf D46): the acquired file, and the derivative the search read when they differ.
                DataFileName = AcquiredFileNameOf(rawFilePath) ?? Path.GetFileName(rawFilePath),
                SearchedDataFileName = AcquiredFileNameOf(rawFilePath) is { } acquired
                    && !string.Equals(acquired, Path.GetFileName(rawFilePath), StringComparison.OrdinalIgnoreCase)
                        ? Path.GetFileName(rawFilePath)
                        : null,
                AssayName = "run " + Path.GetFileNameWithoutExtension(rawFilePath),
                Instrument = ResolveInstrument(rawFilePath),
                PrecursorMassTolerance = common.PrecursorMassTolerance,
                ProductMassTolerance = common.ProductMassTolerance,
                CleavageAgent = common.DigestionParams?.DigestionAgent,
                FixedModifications = ResolveModifications(common.ListOfModsFixed),
                VariableModifications = ResolveModifications(common.ListOfModsVariable),
                DissociationType = common.DissociationType,
                AcquisitionMethod = ResolveAcquisitionMethod(common),
                TechnicalReplicate = technicalReplicate,
                Fraction = fraction
            };

        /// <summary>
        /// One row per annotated channel of an isobaric data file, in reporter m/z order.
        ///
        /// Read from <see cref="TmtFileInfo"/> rather than from mzLib's projection of it, because
        /// <see cref="TmtExperimentalDesign.ToMzLibDesign"/> drops the sample name, and the sample name
        /// is exactly what SDRF's <c>source name</c> is. TmtDesign.txt's replicate and fraction
        /// numbers are already 1-based, so unlike ExperimentalDesign.tsv nothing is added to them.
        ///
        /// EVERY annotated channel is written, including one marked Empty. Skipping Empty channels was
        /// the first attempt and it is wrong twice over (QuantProject thread 002 §5):
        ///
        /// 1. It breaks the round trip BY CONSTRUCTION. The document then has fewer rows than the
        ///    design has channels, so a reader projecting the SDRF back into a design cannot
        ///    reproduce what was searched -- a loss that would have to be enumerated rather than one
        ///    anybody wants.
        /// 2. It changes the plex downstream tooling infers. bigbio's <c>infer_tmtplex</c> reads the
        ///    SET OF LABELS per file; a short set has already been measured inferring the wrong plex,
        ///    and wrong default TMT modifications follow from that.
        ///
        /// An empty channel is a real channel of a real plex, not absent data, and <c>empty</c> is a
        /// value the curated corpus writes. The design does state each channel's sample type
        /// (<see cref="TmtPlexAnnotation.SampleType"/>); writing it is DEFERRED, not impossible.
        /// <c>characteristics[sample type]</c> waits for QuantProject's M4 to put the concept in mzLib
        /// (MAP-08), so every writer spells it one way. When it lands, the only change here is a column
        /// appearing, never rows reappearing or source names changing.
        ///
        /// A channel the design does not annotate at all is still not written, and that is a
        /// different question -- it has no biological replicate, and inventing one is QP-S6, which
        /// QuantProject owes us. <see cref="TmtExperimentalDesign.ToMzLibDesign"/> makes such a
        /// channel Empty for quantification's positional array; doing the same here would put a row
        /// in the file asserting a replicate nobody stated.
        /// </summary>
        private static IEnumerable<SdrfRowInput> BuildChannelRows(TmtFileInfo tmtFile, IsobaricMassTagType tagType,
            CvParam organism, SdrfAssay assay, ISet<string> reusedSampleNames)
        {
            var channelOrder = IsobaricMassTag.GetReporterIonLabels(tagType) ?? new List<string>();
            int OrderOf(string tag)
            {
                int index = channelOrder.FindIndex(l => string.Equals(l, tag?.Trim(), StringComparison.OrdinalIgnoreCase));
                return index < 0 ? int.MaxValue : index;
            }

            string stem = Path.GetFileNameWithoutExtension(tmtFile.FullFilePathWithExtension);
            string plexScope = string.IsNullOrWhiteSpace(tmtFile.Plex) ? stem : tmtFile.Plex.Trim();
            var repeatedInThisPlex = RepeatedWithin(tmtFile.Annotations).ToHashSet(StringComparer.Ordinal);

            foreach (var annotation in tmtFile.Annotations.OrderBy(a => OrderOf(a.Tag)))
            {
                var sample = new SdrfSample
                {
                    // Written, as `not available` when unknown, on channel rows as on file rows (#2816 review).
                    Characteristics = new Dictionary<string, CvParam>
                    {
                        ["characteristics[disease]"] = null,
                        ["characteristics[cell type]"] = null
                    },
                    // The design may leave a channel's sample name blank (an empty channel, or one the
                    // user never named: AnnotatePlexWindow defaults it to ""), and source name is the
                    // column REQ-2 keys sample blocks on -- so it has to be present and distinct. Plex
                    // plus tag is both, and it is the same in every fraction of the plex, which shares
                    // one channel-to-sample map: keyed on the file, one channel of a 12-fraction plex
                    // became 12 samples. The file stem stands in only when the plex has no name. It
                    // deliberately does not spell "empty": that belongs in characteristics[sample type]
                    // when M4 lands, and encoding it here would mean rewriting source names -- and
                    // breaking joins -- on the day it does.
                    //
                    // A name the design reuses in another plex for a different sample is scoped the same way
                    // (SampleNamesReusedForDifferentSamples); a bridge keeps its shared name. A name this
                    // plex gives to several channels for different samples also takes the channel, since
                    // the plex alone would still be one source name (SampleNamesRepeatedWithinAPlex).
                    SourceName = string.IsNullOrWhiteSpace(annotation.SampleName)
                        ? plexScope + " " + annotation.Tag.Trim()
                        : repeatedInThisPlex.Contains(annotation.SampleName)
                            ? plexScope + " " + annotation.SampleName + " " + annotation.Tag.Trim()
                            : reusedSampleNames?.Contains(annotation.SampleName) == true
                                ? plexScope + " " + annotation.SampleName
                                : annotation.SampleName,
                    Organism = organism,
                    BiologicalReplicate = annotation.BiologicalReplicate,
                    Label = ResolveChannelLabel(tagType, annotation.Tag),
                    FactorValue = annotation.Condition,
                    FactorValueColumn = string.IsNullOrWhiteSpace(annotation.Condition)
                        ? null
                        : "factor value[condition]"
                };

                yield return new SdrfRowInput(sample, assay);
            }
        }

        /// <summary>
        /// The parsed TmtDesign.txt keyed by full file path, or null when this is not an isobaric
        /// search or its design cannot be used. Null makes the caller fall back to one row per file.
        ///
        /// Read again rather than shared with multiplex quantification, which keeps its parse in a
        /// local: the file is small, and reading it here leaves the quantification path untouched.
        /// No warning is raised on failure -- quantification has already named the problem, and
        /// SearchTask.WarnAboutSdrfGaps named its consequence for the SDRF before the run.
        /// </summary>
        private Dictionary<string, TmtFileInfo> ReadIsobaricDesignIfPresent(out IsobaricMassTagType? tagType)
        {
            tagType = null;
            if (!Parameters.SearchParameters.DoMultiplexQuantification || Parameters.CurrentRawFileList.Count == 0)
                return null;

            tagType = IsobaricMassTag.GetTagTypeFromModificationId(Parameters.SearchParameters.MultiplexModId);
            if (tagType is null)
                return null;

            string designPath = Path.Combine(
                Path.GetDirectoryName(Parameters.CurrentRawFileList.First()) ?? "",
                GlobalVariables.TmtExperimentalDesignFileName);
            if (!File.Exists(designPath))
                return null;

            var files = TmtExperimentalDesign.Read(designPath, Parameters.CurrentRawFileList, out var errors);
            if (errors.Any())
                return null;

            var byPath = new Dictionary<string, TmtFileInfo>(StringComparer.OrdinalIgnoreCase);
            foreach (var file in files)
                byPath[Path.GetFullPath(file.FullFilePathWithExtension)] = file;
            return byPath;
        }

        /// <summary>
        /// A reporter channel as a PRIDE term: TMT11's "127N" is <c>TMT127N</c>, PRIDE:0000519.
        ///
        /// The design file writes the channel bare, so only the reagent family is added, and it comes
        /// from the search's own tag type. The term itself is looked up in mzLib's pinned PRIDE
        /// vocabulary, never constructed, so a channel PRIDE does not define resolves to nothing. That
        /// is every DiLeu channel, which PRIDE has no terms for.
        ///
        /// Most channels are one PRIDE term whatever the kit: TMT6, TMT10, TMT11 and TMTpro all write
        /// 126 as TMT126. The 131 channel is the exception. PRIDE's TMT131 (PRIDE:0000529) is part of
        /// TMT6PLEX, TMT10PLEX and TMT11PLEX, while TMT131N (PRIDE:0000580) is part of the TMTpro plexes
        /// only. mzLib names TMT10's last channel "131N", so written as given a TMT10 file would carry
        /// a TMTpro term, and tools that infer the plex from the set of labels (sdrf-pipelines'
        /// _infer_plex) would read it as tmt11plex. TMT10's 131N is therefore written TMT131.
        /// </summary>
        private static CvParam ResolveChannelLabel(IsobaricMassTagType tagType, string channel)
        {
            string family = PrideLabelFamily(tagType);

            if (family is null || string.IsNullOrWhiteSpace(channel))
                return null;

            string name = family + channel.Trim();
            if (tagType == IsobaricMassTagType.TMT10 && string.Equals(channel.Trim(), "131N", StringComparison.OrdinalIgnoreCase))
                name = "TMT131";
            // TMT11 keeps TMT131N and TMT131C. PRIDE would put TMT131 in TMT11PLEX too, but sdrf-pipelines'
            // channel_map.yaml spells tmt11plex with 131N/131C, and a 131/131C file matches none of its tables.

            return ControlledVocabulary.Pride.TryGetByName(name, out var term) ? term : null;
        }

        /// <summary>The prefix PRIDE's channel terms carry for this tag type, or null when PRIDE has none.</summary>
        private static string PrideLabelFamily(IsobaricMassTagType tagType) => tagType switch
        {
            IsobaricMassTagType.TMT6 or IsobaricMassTagType.TMT10 or IsobaricMassTagType.TMT11
                or IsobaricMassTagType.TMT18 => "TMT",
            IsobaricMassTagType.iTRAQ4 or IsobaricMassTagType.iTRAQ8 => "ITRAQ",
            _ => null
        };

        /// <summary>
        /// The acquired file a searched file came from: the run's starting file with the same stem once the
        /// -calib / -averaged suffixes are removed, or null when there is no starting list (a task run on its own)
        /// or no starting file matches.
        /// </summary>
        private string AcquiredFileNameOf(string searchedPath)
        {
            if (Parameters.AcquiredSpectraFiles is not { Count: > 0 } acquired) return null;
            string stem = Path.GetFileNameWithoutExtension(searchedPath);
            string previous;
            do
            {
                previous = stem;
                foreach (string suffix in new[] { CalibrationTask.CalibSuffix, SpectralAveragingTask.AveragingSuffix })
                    if (stem.EndsWith(suffix, StringComparison.OrdinalIgnoreCase))
                        stem = stem[..^suffix.Length];
            } while (stem != previous);
            string match = acquired.FirstOrDefault(f =>
                string.Equals(Path.GetFileNameWithoutExtension(f), stem, StringComparison.OrdinalIgnoreCase));
            return match is null ? null : Path.GetFileName(match);
        }

        /// <summary>
        /// The experimental design keyed by file stem, or empty when the user did not supply one.
        /// Never the fabricated fallback: see the remarks in <see cref="BuildSdrfRows"/>.
        /// </summary>
        private Dictionary<string, SpectraFileInfo> ReadExperimentalDesignIfPresent()
        {
            var byStem = new Dictionary<string, SpectraFileInfo>(StringComparer.OrdinalIgnoreCase);
            if (Parameters.CurrentRawFileList.Count == 0) return byStem;

            string designPath = Path.Combine(
                Path.GetDirectoryName(Parameters.CurrentRawFileList.First()) ?? "",
                GlobalVariables.ExperimentalDesignFileName);

            if (!File.Exists(designPath))
            {
                Warn($"No {GlobalVariables.ExperimentalDesignFileName} beside the spectra files, so the " +
                     "SDRF cannot describe conditions, replicates or fractions. It will record the " +
                     "search parameters and what can be read from the files themselves; replicate and " +
                     "fraction are written as 'not available'.");
                return byStem;
            }

            List<SpectraFileInfo> infos;
            List<string> errors;
            try
            {
                infos = ExperimentalDesign.ReadExperimentalDesign(designPath, Parameters.CurrentRawFileList, out errors);
            }
            catch (Exception e) when (e is IOException or UnauthorizedAccessException)
            {
                // Open in Excel, say. The pre-run warning promised the SDRF would then describe the search
                // only; failing here would lose the whole SDRF instead.
                Warn($"{GlobalVariables.ExperimentalDesignFileName} could not be read ({e.Message}), so the SDRF " +
                     "describes the search only: no conditions, replicates or fractions.");
                return byStem;
            }
            if (errors.Any())
            {
                Warn($"{GlobalVariables.ExperimentalDesignFileName} has errors, so it is not being used " +
                     $"for the SDRF: {string.Join("; ", errors)}");
                return byStem;
            }

            foreach (var info in infos)
                byStem[info.FilenameWithoutExtension] = info;
            return byStem;
        }

        /// <summary>
        /// The organism, taken from the search database rather than looked up, so no taxonomy
        /// ontology has to be shipped or queried. A UniProt database states it and mzLib retains it
        /// as <c>Protein.NcbiTaxonomyId</c>.
        ///
        /// Today that reaches a MetaMorpheus search only from a UniProt XML database. mzLib can read
        /// FASTA's OX= too, but MetaMorpheus passes LoadProteinFasta explicit regexes without the
        /// organism-id one, so a FASTA search leaves the column "not available". That loader call is
        /// #2782's to change, and it deliberately has not.
        /// </summary>
        private CvParam ResolveOrganismFromSearchDatabase()
        {
            // Never a contaminant's taxon: MetaMorpheusContaminants.xml carries NCBI Taxonomy 9913, and a FASTA
            // target carries none, so the commonest search would otherwise say every sample is Bos taurus. And
            // only when the targets name ONE organism: "first protein wins" would pick one of several at random.
            var taxa = Parameters.BioPolymerList?
                .OfType<Protein>()
                .Where(p => !p.IsDecoy && !p.IsContaminant && !string.IsNullOrEmpty(p.NcbiTaxonomyId))
                .GroupBy(p => p.NcbiTaxonomyId)
                .ToList();

            if (taxa is null || taxa.Count != 1) return null;

            var first = taxa[0].First();
            return new CvParam("NCBITaxon", "NCBITaxon:" + first.NcbiTaxonomyId, first.Organism ?? "", "");
        }

        /// <summary>
        /// The instrument, as the search's own load of the file read it (SourceFile). mzML carries an accessioned
        /// term; a Thermo RAW carries only a name, which SdrfBuilder resolves against PSI-MS. Opening every file
        /// again here would re-read the whole dataset, and a RAW's SourceFile hashes the entire file first, so a
        /// file the search did not record resolves to no instrument (written as `not available`).
        /// </summary>
        private CvParam ResolveInstrument(string rawFilePath) =>
            Parameters.InstrumentModelsByFile is { } models && models.TryGetValue(rawFilePath, out var model) ? model : null;

        /// <summary>
        /// Only a genuinely label-free search is described as label free.
        ///
        /// This resolves the label for a row that describes a whole FILE. SDRF wants one row per sample
        /// per channel for a labelled run. An isobaric run with a usable TmtDesign.txt gets those rows
        /// from <see cref="BuildChannelRows"/> and never reaches here; one without falls back to a row
        /// per file. SILAC has no channel-to-sample mapping in MetaMorpheus at all. Either way the
        /// label is left unresolved and the coverage report shows it; guessing would invent an
        /// experimental design.
        ///
        /// Returning "label free sample" for a labelled run would be worse than returning nothing:
        /// the column comes out fully populated with a confident falsehood, which
        /// <see cref="SdrfCoverage"/> cannot flag because it only measures emptiness.
        ///
        /// Static and parameterised so the decision can be tested without driving a whole search.
        /// </summary>
        private static CvParam ResolveLabel(SearchParameters searchParameters)
        {
            if (searchParameters.DoMultiplexQuantification)
                return null;

            // SILAC, including the turnover variants, which carry their labels separately.
            if (searchParameters.SilacLabels?.Any() == true
                || searchParameters.StartTurnoverLabel is not null
                || searchParameters.EndTurnoverLabel is not null)
                return null;

            return new CvParam("MS", "MS:1002038", "label free sample", "");
        }

        /// <summary>
        /// The acquisition method, read from the search rather than assumed.
        ///
        /// Every MetaMorpheus search today is data-dependent, which is exactly why this was
        /// hardcoded and exactly why that was a hazard: a constant is not wrong until the day the
        /// capability lands, and then it is wrong invisibly. The column would be 100% filled with a
        /// false CV term, so no coverage or drift instrument could see it.
        ///
        /// In-source decay deliberately resolves to nothing. It is not a precursor-selection scheme
        /// and PSI-MS/PRIDE define no acquisition-method term for it; borrowing the nearest-looking
        /// one is the misannotation D12 exists to keep out of authored output.
        /// </summary>
        private static CvParam ResolveAcquisitionMethod(CommonParameters commonParameters)
        {
            if (commonParameters?.DIAparameters is null)
                return new CvParam("PRIDE", "PRIDE:0000627", "Data-dependent acquisition", "");

            return commonParameters.DIAparameters.AanalysisType switch
            {
                DIAanalysisType.DIA => new CvParam("PRIDE", "PRIDE:0000450", "Data-independent acquisition", ""),
                _ => null
            };
        }

        private static IReadOnlyList<Modification> ResolveModifications(IEnumerable<(string, string)> mods)
        {
            if (mods is null) return Array.Empty<Modification>();

            // The (ModificationType, IdWithMotif) tuples are resolved to real Modification objects
            // so the builder can read the UNIMOD accession off DatabaseReference rather than being
            // handed a string to guess from. Both halves are the identity: the same IdWithMotif is
            // known under more than one type (Hex on Y is one), and keying on
            // the id alone would write every one of them.
            var wanted = new HashSet<(string, string)>(mods);
            return GlobalVariables.AllModsKnown
                .Where(m => wanted.Contains((m.ModificationType, m.IdWithMotif)))
                .ToList();
        }

        /// <summary>
        /// Reports which columns the document did not actually populate.
        ///
        /// Without this an SDRF full of "not available" looks like a success: the validator accepts
        /// reserved words as the correct way to state an absence, and the drift lint skips them, so
        /// nothing else in the stack can see an empty corpus.
        /// </summary>
        private void ReportSdrfCoverage(SdrfDocument document)
        {
            var uninformative = SdrfCoverage
                .Uninformative(new SdrfCollection(new[] { document }, new[] { "this run" }))
                .Where(c => c.FillRate < 1.0)
                .Select(c => c.Column)
                .ToList();

            string designFileName = Parameters.SearchParameters.DoMultiplexQuantification
                ? GlobalVariables.TmtExperimentalDesignFileName
                : GlobalVariables.ExperimentalDesignFileName;

            if (uninformative.Any())
                Warn("The SDRF was written, but these columns say nothing and cannot be mined: " +
                     string.Join(", ", uninformative) + ". " + designFileName +
                     " can supply conditions, biological and technical replicates and fractions; everything " +
                     "else (organism part, disease, cell type and the like) has to be curated by hand in the " +
                     "written SDRF, since MetaMorpheus reads no input SDRF.");
        }
    }
}
