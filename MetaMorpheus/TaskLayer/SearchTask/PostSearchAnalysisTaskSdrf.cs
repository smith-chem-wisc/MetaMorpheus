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
            // which every file is its own biological replicate with an empty condition. Emitting
            // that as though it were curated metadata is the worst thing this writer could do -- it
            // looks like data. So read the file, and if it is not there say so and describe only
            // what is actually known.
            //
            // An isobaric search keeps its design in TmtDesign.txt instead, and that file is the only
            // place the channel-to-sample map exists. When it is usable the SDRF gets one row per
            // sample per channel, as the specification wants; when it is not, the search is described
            // one row per file with the label left unresolved, exactly as before.
            var isobaric = ReadIsobaricDesignIfPresent(out var tagType);
            var design = isobaric is null
                ? ReadExperimentalDesignIfPresent()
                : new Dictionary<string, SpectraFileInfo>(StringComparer.OrdinalIgnoreCase);

            var organism = ResolveOrganismFromSearchDatabase();

            foreach (var rawFilePath in Parameters.CurrentRawFileList)
            {
                // The ORIGINAL file, not a -calib or -averaged derivative: those are intermediates
                // this run produced, and an SDRF describes the data as acquired.
                string fileName = Path.GetFileName(rawFilePath);
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

                if (isobaric is not null && isobaric.TryGetValue(Path.GetFullPath(rawFilePath), out var tmtFile))
                {
                    foreach (var row in BuildChannelRows(tmtFile, tagType!.Value, organism,
                                 BuildAssay(rawFilePath, common, tmtFile.TechnicalReplicate, tmtFile.Fraction)))
                        yield return row;
                    continue;
                }

                design.TryGetValue(stem, out var sampleInfo);

                var sample = new SdrfSample
                {
                    SourceName = sampleInfo?.Condition is { Length: > 0 } condition
                        ? condition + " " + (sampleInfo.BiologicalReplicate + 1)
                        : stem,
                    Organism = organism,
                    // SDRF is 1-based; SpectraFileInfo stores these 0-based.
                    BiologicalReplicate = (sampleInfo?.BiologicalReplicate ?? 0) + 1,
                    Label = ResolveLabel(Parameters.SearchParameters),
                    FactorValue = sampleInfo?.Condition,
                    FactorValueColumn = string.IsNullOrWhiteSpace(sampleInfo?.Condition)
                        ? null
                        : "factor value[condition]"
                };

                var assay = BuildAssay(rawFilePath, common,
                    technicalReplicate: (sampleInfo?.TechnicalReplicate ?? 0) + 1,
                    fraction: (sampleInfo?.Fraction ?? 0) + 1);

                yield return new SdrfRowInput(sample, assay);
            }
        }

        /// <summary>
        /// The assay half of a row: everything the search itself knows about one data file. Shared by
        /// every channel row of an isobaric file, which is what makes those rows one assay.
        /// </summary>
        private SdrfAssay BuildAssay(string rawFilePath, CommonParameters common, int technicalReplicate, int fraction) =>
            new SdrfAssay
            {
                DataFileName = Path.GetFileName(rawFilePath),
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
        /// value the curated corpus writes. Still no <c>characteristics[sample type]</c> column until
        /// QuantProject's M4 puts the concept in mzLib: when it lands, the only change here is a
        /// column appearing, never rows reappearing.
        ///
        /// A channel the design does not annotate at all is still not written, and that is a
        /// different question -- it has no biological replicate, and inventing one is QP-S6, which
        /// QuantProject owes us. <see cref="TmtExperimentalDesign.ToMzLibDesign"/> makes such a
        /// channel Empty for quantification's positional array; doing the same here would put a row
        /// in the file asserting a replicate nobody stated.
        /// </summary>
        private static IEnumerable<SdrfRowInput> BuildChannelRows(TmtFileInfo tmtFile, IsobaricMassTagType tagType,
            CvParam organism, SdrfAssay assay)
        {
            var channelOrder = IsobaricMassTag.GetReporterIonLabels(tagType) ?? new List<string>();
            int OrderOf(string tag)
            {
                int index = channelOrder.FindIndex(l => string.Equals(l, tag?.Trim(), StringComparison.OrdinalIgnoreCase));
                return index < 0 ? int.MaxValue : index;
            }

            string stem = Path.GetFileNameWithoutExtension(tmtFile.FullFilePathWithExtension);

            foreach (var annotation in tmtFile.Annotations.OrderBy(a => OrderOf(a.Tag)))
            {
                var sample = new SdrfSample
                {
                    // An empty channel is the one case where the design may leave the sample name
                    // blank, and source name is the column REQ-2 keys sample blocks on -- so it has
                    // to be present and distinct. File stem plus tag is both, and it is built only
                    // from facts already in the row. It deliberately does not spell "empty": that
                    // belongs in characteristics[sample type] when M4 lands, and encoding it here
                    // would mean rewriting source names -- and breaking joins -- on the day it does.
                    SourceName = string.IsNullOrWhiteSpace(annotation.SampleName)
                        ? stem + " " + annotation.Tag.Trim()
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
                Directory.GetParent(Parameters.CurrentRawFileList.First())!.ToString(),
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
        /// TMTpro 18-plex channels share their names with TMT's in PRIDE (TMT126 ... TMT135N); PRIDE
        /// distinguishes the kit, not the channel.
        /// </summary>
        private static CvParam ResolveChannelLabel(IsobaricMassTagType tagType, string channel)
        {
            string family = tagType switch
            {
                IsobaricMassTagType.TMT6 or IsobaricMassTagType.TMT10 or IsobaricMassTagType.TMT11
                    or IsobaricMassTagType.TMT18 => "TMT",
                IsobaricMassTagType.iTRAQ4 or IsobaricMassTagType.iTRAQ8 => "ITRAQ",
                _ => null
            };

            if (family is null || string.IsNullOrWhiteSpace(channel))
                return null;

            return ControlledVocabulary.Pride.TryGetByName(family + channel.Trim(), out var term) ? term : null;
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
                Directory.GetParent(Parameters.CurrentRawFileList.First())!.ToString(),
                GlobalVariables.ExperimentalDesignFileName);

            if (!File.Exists(designPath))
            {
                Warn($"No {GlobalVariables.ExperimentalDesignFileName} beside the spectra files, so the " +
                     "SDRF cannot describe conditions, replicates or fractions. It will record the " +
                     "search parameters and what can be read from the files themselves.");
                return byStem;
            }

            var infos = ExperimentalDesign.ReadExperimentalDesign(
                designPath, Parameters.CurrentRawFileList, out var errors);
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
            var withTaxon = Parameters.BioPolymerList?
                .OfType<Protein>()
                .FirstOrDefault(p => !p.IsDecoy && !string.IsNullOrEmpty(p.NcbiTaxonomyId));

            if (withTaxon is null) return null;

            return new CvParam("NCBITaxon", "NCBITaxon:" + withTaxon.NcbiTaxonomyId,
                withTaxon.Organism ?? "", "");
        }

        /// <summary>
        /// The instrument, read from the data file itself. mzML carries an accessioned term; a
        /// Thermo RAW carries only a name, which SdrfBuilder resolves against PSI-MS.
        /// </summary>
        private CvParam ResolveInstrument(string rawFilePath)
        {
            try
            {
                return MsDataFileReader.GetDataFile(rawFilePath).GetSourceFile()?.InstrumentModel;
            }
            catch (Exception)
            {
                // Reading a header must never be the thing that fails a completed search.
                return null;
            }
        }

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

            return new CvParam("", "", "label free sample", "");
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
                     string.Join(", ", uninformative) + ". Supply them in " +
                     designFileName + " or an input SDRF.");
        }
    }
}
