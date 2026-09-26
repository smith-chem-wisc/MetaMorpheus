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
                    RequireSampleMetadata = false
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
            // not read here. Without a design the replicate and fraction numbers are UNKNOWN, but
            // SdrfSample / SdrfAssay in the mzLib this builds against take a plain int, so they are
            // written as 1. mzLib #1378 makes them nullable ("not available"); pass null here once a
            // release carries it. Until then a design-less SDRF states 1 for each, and the warning
            // below says so.
            var design = ReadExperimentalDesignIfPresent();

            var organism = ResolveOrganismFromSearchDatabase();

            foreach (var rawFilePath in Parameters.CurrentRawFileList)
            {
                // Two names (sdrf D46). comment[data file] is the ACQUIRED file, as deposited -- what the
                // SDRF specification means by the column and what SdrfAssay.DataFileName documents. After a
                // Calibrate or Average task the search reads a -calib / -averaged derivative instead, and that
                // name goes to comment[searched data file]. The acquired name, extension included, comes from
                // the files the run started from; a task run on its own was given the acquired files.
                string searchedName = Path.GetFileName(rawFilePath);
                string acquiredName = AcquiredFileNameOf(rawFilePath) ?? searchedName;
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
                    // SDRF is 1-based; SpectraFileInfo stores these 0-based.
                    BiologicalReplicate = (sampleInfo?.BiologicalReplicate ?? 0) + 1,
                    Label = ResolveLabel(Parameters.SearchParameters),
                    FactorValue = sampleInfo?.Condition,
                    FactorValueColumn = string.IsNullOrWhiteSpace(sampleInfo?.Condition)
                        ? null
                        : "factor value[condition]"
                };

                var assay = new SdrfAssay
                {
                    DataFileName = acquiredName,
                    SearchedDataFileName = string.Equals(acquiredName, searchedName, StringComparison.OrdinalIgnoreCase)
                        ? null
                        : searchedName,
                    AssayName = "run " + stem,
                    Instrument = ResolveInstrument(rawFilePath),
                    PrecursorMassTolerance = common.PrecursorMassTolerance,
                    ProductMassTolerance = common.ProductMassTolerance,
                    CleavageAgent = common.DigestionParams?.DigestionAgent,
                    FixedModifications = ResolveModifications(common.ListOfModsFixed),
                    VariableModifications = ResolveModifications(common.ListOfModsVariable),
                    DissociationType = common.DissociationType,
                    AcquisitionMethod = ResolveAcquisitionMethod(common),
                    TechnicalReplicate = (sampleInfo?.TechnicalReplicate ?? 0) + 1,
                    Fraction = (sampleInfo?.Fraction ?? 0) + 1
                };

                yield return new SdrfRowInput(sample, assay);
            }
        }

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
                     "search parameters and what can be read from the files themselves, and it writes " +
                     "replicate and fraction 1 for every file, which nothing established.");
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
        /// SDRF wants one row per sample per CHANNEL, and this writer emits one row per file. SILAC
        /// has no channel-to-sample mapping in MetaMorpheus at all. Isobaric runs do, in
        /// TmtDesign.txt, but expanding rows from it is not done here yet. Until it is, the label is
        /// left unresolved and the coverage report shows it; guessing would invent an experimental
        /// design.
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

            if (uninformative.Any())
                Warn("The SDRF was written, but these columns say nothing and cannot be mined: " +
                     string.Join(", ", uninformative) + ". " + GlobalVariables.ExperimentalDesignFileName +
                     " can supply conditions, biological and technical replicates and fractions; everything " +
                     "else (organism part, disease, cell type and the like) has to be curated by hand in the " +
                     "written SDRF, since MetaMorpheus reads no input SDRF.");
        }
    }
}
