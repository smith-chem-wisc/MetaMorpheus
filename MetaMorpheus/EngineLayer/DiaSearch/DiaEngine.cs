// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using MassSpectrometry;
using MassSpectrometry.Dia;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.DiaSearch
{
    /// <summary>
    /// MetaMorpheus engine facade for DIA spectral library search.
    /// 
    /// Orchestrates the full mzLib DIA pipeline:
    ///   1. Convert LibrarySpectrum → LibraryPrecursorInput (bridge Omics → mzLib)
    ///   2. Build DiaScanIndex from raw scans (SoA layout)
    ///   3. Generate fragment queries from library precursors
    ///   4. Create fragment extractor (CPU or GPU)
    ///   5. Extract all fragment XICs in parallel
    ///   6. Assemble and score results
    /// 
    /// Implements IDisposable to clean up DiaScanIndex and DiaExtractionOrchestrator.
    /// </summary>
    public class DiaEngine : MetaMorpheusEngine, IDisposable
    {
        private readonly MsDataScan[] _scans;
        private readonly List<LibrarySpectrum> _librarySpectra;
        private readonly DiaSearchParameters _diaParams;
        private DiaScanIndex _scanIndex;
        private DiaExtractionOrchestrator _orchestrator;
        private bool _disposed;

        /// <summary>
        /// Creates a new DIA search engine.
        /// </summary>
        /// <param name="scans">All scans from the raw file (MS1 + MS2; MS2 are filtered during indexing).</param>
        /// <param name="librarySpectra">Spectral library entries to search against.</param>
        /// <param name="diaParams">mzLib-level DIA search parameters.</param>
        /// <param name="commonParameters">Standard MetaMorpheus common parameters.</param>
        /// <param name="fileSpecificParameters">File-specific parameter overrides.</param>
        /// <param name="nestedIds">Nested task identifiers for status reporting.</param>
        public DiaEngine(
            MsDataScan[] scans,
            List<LibrarySpectrum> librarySpectra,
            DiaSearchParameters diaParams,
            CommonParameters commonParameters,
            List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
            List<string> nestedIds)
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            _scans = scans ?? throw new ArgumentNullException(nameof(scans));
            _librarySpectra = librarySpectra ?? throw new ArgumentNullException(nameof(librarySpectra));
            _diaParams = diaParams ?? throw new ArgumentNullException(nameof(diaParams));
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            // ── Step 1: Convert LibrarySpectrum → LibraryPrecursorInput ────
            Status("Converting library spectra...");
            var precursorInputs = ConvertLibrarySpectra(_librarySpectra);

            if (precursorInputs.Length == 0)
            {
                return new DiaEngineResults(this, new List<DiaSearchResult>());
            }

            // ── Step 2: Build SoA scan index ───────────────────────────────
            Status("Building DIA scan index...");
            _scanIndex = DiaScanIndexBuilder.Build(_scans);

            if (_scanIndex.ScanCount == 0)
            {
                return new DiaEngineResults(this, new List<DiaSearchResult>());
            }

            // ── Step 3: Generate fragment queries ──────────────────────────
            Status("Generating fragment queries...");
            var generationResult = DiaLibraryQueryGenerator.Generate(precursorInputs, _scanIndex, _diaParams);

            if (generationResult.Queries.Length == 0)
            {
                return new DiaEngineResults(this, new List<DiaSearchResult>());
            }

            // ── Step 4: Create fragment extractor factory ──────────────────
            var extractorFactory = FragmentExtractorFactory.CreateFactory(
                preferCpu: !_diaParams.PreferGpu);

            Status($"Extraction backend: {FragmentExtractorFactory.DescribeBackend(!_diaParams.PreferGpu)}");

            // ── Step 5: Extract all fragment XICs ──────────────────────────
            Status($"Extracting {generationResult.Queries.Length:N0} fragment XICs across " +
                   $"{_scanIndex.WindowCount} windows...");

            _orchestrator = new DiaExtractionOrchestrator(_scanIndex, extractorFactory);
            var extractionResult = _orchestrator.ExtractAll(
                generationResult.Queries,
                _diaParams.EffectiveMaxThreads);

            // ── Step 6: Assemble and score results ─────────────────────────
            Status("Scoring and assembling results...");

            var dotProductScorer = new NormalizedDotProductScorer();
            var spectralAngleScorer = new SpectralAngleScorer();

            var diaResults = DiaLibraryQueryGenerator.AssembleResults(
                precursorInputs,
                generationResult,
                extractionResult.Results,
                _diaParams,
                dotProductScorer,
                spectralAngleScorer);

            Status($"DIA search complete: {diaResults.Count(r => !r.IsDecoy):N0} targets, " +
                   $"{diaResults.Count(r => r.IsDecoy):N0} decoys " +
                   $"(skipped: {generationResult.SkippedNoWindow} no window, " +
                   $"{generationResult.SkippedNoFragments} no fragments)");

            return new DiaEngineResults(this, diaResults);
        }

        /// <summary>
        /// Converts LibrarySpectrum objects (Omics layer) into LibraryPrecursorInput structs
        /// (mzLib DIA engine layer). This bridges the Omics dependency boundary so that
        /// the mzLib DIA engine has no direct dependency on Omics types.
        /// </summary>
        /// <param name="librarySpectra">Library spectra from .msp reader or Koina predictions.</param>
        /// <returns>Array of lightweight precursor inputs for the DIA query generator.</returns>
        public static LibraryPrecursorInput[] ConvertLibrarySpectra(List<LibrarySpectrum> librarySpectra)
        {
            if (librarySpectra == null || librarySpectra.Count == 0)
                return Array.Empty<LibraryPrecursorInput>();

            var inputs = new LibraryPrecursorInput[librarySpectra.Count];

            for (int i = 0; i < librarySpectra.Count; i++)
            {
                var lib = librarySpectra[i];
                var fragments = lib.MatchedFragmentIons;

                float[] fragMzs = new float[fragments.Count];
                float[] fragIntensities = new float[fragments.Count];

                for (int f = 0; f < fragments.Count; f++)
                {
                    fragMzs[f] = (float)fragments[f].Mz;
                    fragIntensities[f] = (float)fragments[f].Intensity;
                }

                inputs[i] = new LibraryPrecursorInput(
                    sequence: lib.Sequence,
                    precursorMz: lib.PrecursorMz,
                    chargeState: lib.ChargeState,
                    retentionTime: lib.RetentionTime,
                    isDecoy: lib.IsDecoy,
                    fragmentMzs: fragMzs,
                    fragmentIntensities: fragIntensities);
            }

            return inputs;
        }

        public void Dispose()
        {
            if (_disposed) return;
            _disposed = true;

            _orchestrator?.Dispose();
            _scanIndex?.Dispose();
        }
    }
}
