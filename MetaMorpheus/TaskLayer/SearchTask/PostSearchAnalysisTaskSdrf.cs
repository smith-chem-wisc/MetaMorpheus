using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
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
            // annotation must not destroy a search that has already succeeded. The requirement that
            // the metadata EXISTS is enforced earlier, at task validation, where failing is cheap.
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
                    // The sample metadata is only as complete as the input allowed. Enforcement
                    // happened at validation time; by here the run is committed, so accept what we
                    // have and let SdrfCoverage report on it rather than throwing mid-write.
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
            // which every file is its own biological replicate with an empty condition. Emitting
            // that as though it were curated metadata is the worst thing this writer could do -- it
            // looks like data. So read the file, and if it is not there say so and describe only
            // what is actually known.
            var design = ReadExperimentalDesignIfPresent();

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
                var common = FileSpecificParameters
                    ?.FirstOrDefault(f => string.Equals(f.FileName, fileName, StringComparison.OrdinalIgnoreCase))
                    .Parameters ?? CommonParameters;

                design.TryGetValue(stem, out var sampleInfo);

                var sample = new SdrfSample
                {
                    SourceName = sampleInfo?.Condition is { Length: > 0 } condition
                        ? condition + " " + (sampleInfo.BiologicalReplicate + 1)
                        : stem,
                    Organism = organism,
                    // SDRF is 1-based; SpectraFileInfo stores these 0-based.
                    BiologicalReplicate = (sampleInfo?.BiologicalReplicate ?? 0) + 1,
                    Label = ResolveLabel(),
                    FactorValue = sampleInfo?.Condition,
                    FactorValueColumn = string.IsNullOrWhiteSpace(sampleInfo?.Condition)
                        ? null
                        : "factor value[condition]"
                };

                var assay = new SdrfAssay
                {
                    DataFileName = fileName,
                    AssayName = "run " + stem,
                    Instrument = ResolveInstrument(rawFilePath),
                    PrecursorMassTolerance = common.PrecursorMassTolerance,
                    ProductMassTolerance = common.ProductMassTolerance,
                    CleavageAgent = common.DigestionParams?.DigestionAgent,
                    FixedModifications = ResolveModifications(common.ListOfModsFixed),
                    VariableModifications = ResolveModifications(common.ListOfModsVariable),
                    DissociationType = common.DissociationType,
                    AcquisitionMethod = new CvParam("PRIDE", "PRIDE:0000627", "Data-dependent acquisition", ""),
                    TechnicalReplicate = (sampleInfo?.TechnicalReplicate ?? 0) + 1,
                    Fraction = (sampleInfo?.Fraction ?? 0) + 1
                };

                yield return new SdrfRowInput(sample, assay);
            }
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
        /// The organism, taken from the search database rather than looked up. A UniProt database
        /// states it, and mzLib retains it as <c>Protein.NcbiTaxonomyId</c> for both the XML and
        /// FASTA paths, so no taxonomy ontology has to be shipped or queried.
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
        /// Label-free and SILAC are expressible; isobaric labelling is not, and says so.
        ///
        /// SDRF wants one row per sample per CHANNEL, and MetaMorpheus models a single tag type for
        /// the whole search with no channel-to-sample mapping anywhere. Guessing one would invent
        /// an experimental design, so the label is left unresolved and the coverage report shows it.
        /// </summary>
        private CvParam ResolveLabel()
        {
            if (Parameters.SearchParameters.DoMultiplexQuantification)
                return null;

            return new CvParam("", "", "label free sample", "");
        }

        private IReadOnlyList<Modification> ResolveModifications(IEnumerable<(string, string)> mods)
        {
            if (mods is null) return Array.Empty<Modification>();

            // The (ModificationType, IdWithMotif) tuples are resolved to real Modification objects
            // so the builder can read the UNIMOD accession off DatabaseReference rather than being
            // handed a string to guess from.
            var wanted = new HashSet<string>(mods.Select(m => m.Item2), StringComparer.Ordinal);
            return GlobalVariables.AllModsKnown
                .Where(m => wanted.Contains(m.IdWithMotif))
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
                     string.Join(", ", uninformative) + ". Supply them in " +
                     GlobalVariables.ExperimentalDesignFileName + " or an input SDRF.");
        }
    }
}
