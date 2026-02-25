// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Dia;
using Omics.SpectrumMatch;

namespace EngineLayer.DiaSearch
{
    /// <summary>
    /// MetaMorpheus-side DIA search engine.
    /// 
    /// Extends MetaMorpheusEngine so it integrates with the standard Run()/RunSpecific()
    /// pattern, event handlers, timing, and task orchestration.
    /// 
    /// Bridges the Omics layer (LibrarySpectrum, MatchedFragmentIon) with the mzLib DIA
    /// extraction engine (DiaScanIndex, DiaLibraryQueryGenerator, DiaExtractionOrchestrator,
    /// Scorers). The bridge is the ConvertLibrarySpectra() method which converts
    /// LibrarySpectrum → LibraryPrecursorInput (dependency-free struct in mzLib).
    /// 
    /// Pipeline (executed in RunSpecific):
    ///   1. Convert LibrarySpectrum list → LibraryPrecursorInput[]
    ///   2. Build DiaScanIndex from MsDataScan[]
    ///   3. Generate FragmentQuery[] via DiaLibraryQueryGenerator
    ///   4. Extract via DiaExtractionOrchestrator (parallel, CPU or GPU)
    ///   5. Score and assemble DiaSearchResult list
    ///   6. Return DiaEngineResults
    /// </summary>
    public sealed class DiaEngine : MetaMorpheusEngine, IDisposable
    {
        private readonly MsDataScan[] _scans;
        private readonly List<LibrarySpectrum> _librarySpectra;
        private readonly DiaSearchParameters _diaParameters;
        private DiaScanIndex _index;
        private bool _disposed;

        /// <summary>
        /// Creates a DIA search engine following the MetaMorpheusEngine pattern.
        /// </summary>
        /// <param name="scans">All scans from the DIA data file (MS1 + MS2; MS1 filtered internally).</param>
        /// <param name="librarySpectra">Spectral library entries to search against.</param>
        /// <param name="diaParameters">DIA-specific search parameters.</param>
        /// <param name="commonParameters">Standard MetaMorpheus common parameters.</param>
        /// <param name="fileSpecificParameters">Per-file parameter overrides.</param>
        /// <param name="nestedIds">Task/engine identification for event reporting.</param>
        public DiaEngine(
            MsDataScan[] scans,
            List<LibrarySpectrum> librarySpectra,
            DiaSearchParameters diaParameters,
            CommonParameters commonParameters,
            List<(string FileName, CommonParameters Parameters)> fileSpecificParameters,
            List<string> nestedIds)
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            _scans = scans ?? throw new ArgumentNullException(nameof(scans));
            _librarySpectra = librarySpectra ?? throw new ArgumentNullException(nameof(librarySpectra));
            _diaParameters = diaParameters ?? new DiaSearchParameters();
        }

        /// <summary>
        /// Executes the DIA search pipeline.
        /// Called by MetaMorpheusEngine.Run(), which handles timing and events.
        /// </summary>
        protected override MetaMorpheusEngineResults RunSpecific()
        {
            // ── Convert library ─────────────────────────────────────────────
            Status("Converting spectral library...");
            var precursorInputs = ConvertLibrarySpectra(_librarySpectra);

            int targetCount = 0;
            int decoyCount = 0;
            for (int i = 0; i < precursorInputs.Length; i++)
            {
                if (precursorInputs[i].IsDecoy)
                    decoyCount++;
                else
                    targetCount++;
            }

            // Early exit: empty library
            if (precursorInputs.Length == 0)
            {
                return new DiaEngineResults(this,
                    new List<DiaSearchResult>(),
                    targetLibraryCount: 0,
                    decoyLibraryCount: 0);
            }

            // ── Step 1: Build SoA scan index ────────────────────────────────
            Status("Building DIA scan index...");
            _index = DiaScanIndexBuilder.Build(_scans);

            Status($"Indexed {_index.ScanCount} MS2 scans across {_index.WindowCount} windows " +
                   $"({_index.TotalPeakCount:N0} peaks)");

            // ── Step 2: Generate fragment queries ───────────────────────────
            Status("Generating fragment queries...");
            var generationResult = DiaLibraryQueryGenerator.Generate(
                precursorInputs, _index, _diaParameters);

            Status($"Generated {generationResult.Queries.Length:N0} queries " +
                   $"({generationResult.SkippedNoWindow} precursors outside windows)");

            // Early exit: no valid queries (all precursors outside windows)
            if (generationResult.Queries.Length == 0)
            {
                _index.Dispose();
                _index = null;
                return new DiaEngineResults(this,
                    new List<DiaSearchResult>(),
                    targetLibraryCount: targetCount,
                    decoyLibraryCount: decoyCount);
            }

            // ── Step 3: Extract fragment XICs (parallel) ────────────────────
            Status("Extracting fragment ion chromatograms...");
            var extractorFactory = FragmentExtractorFactory.CreateFactory(
                preferCpu: !_diaParameters.PreferGpu);

            using var orchestrator = new DiaExtractionOrchestrator(_index, extractorFactory);
            var extractionResult = orchestrator.ExtractAll(
                generationResult.Queries,
                maxDegreeOfParallelism: _diaParameters.EffectiveMaxThreads);

            Status($"Extracted {extractionResult.TotalDataPoints:N0} XIC data points");

            // ── Step 4: Score and assemble results ──────────────────────────
            Status("Scoring precursor matches...");
            var dotProductScorer = new NormalizedDotProductScorer();
            var spectralAngleScorer = new SpectralAngleScorer();

            var results = DiaLibraryQueryGenerator.AssembleResults(
                precursorInputs,
                generationResult,
                extractionResult.Results,
                _diaParameters,
                dotProductScorer,
                spectralAngleScorer);

            // Report summary
            int targetHits = 0;
            int decoyHits = 0;
            for (int i = 0; i < results.Count; i++)
            {
                if (results[i].IsDecoy)
                    decoyHits++;
                else
                    targetHits++;
            }

            Status($"DIA search complete: {results.Count:N0} precursor matches " +
                   $"({targetHits} targets, {decoyHits} decoys)");

            // Clean up index
            _index.Dispose();
            _index = null;

            return new DiaEngineResults(this, results,
                targetLibraryCount: targetCount,
                decoyLibraryCount: decoyCount);
        }

        /// <summary>
        /// Converts Omics LibrarySpectrum objects into mzLib-native LibraryPrecursorInput
        /// structs. This is the bridge between MetaMorpheus (which has access to Omics types)
        /// and the mzLib DIA engine (which is Omics-independent).
        /// 
        /// Key conversions:
        ///   - PrecursorMz: kept as double (precision matters for window matching)
        ///   - RetentionTime: kept as double? (precision matters for RT windowing)
        ///   - FragmentMzs/Intensities: converted to float[] (sufficient for MS2 fragments)
        ///   - Skips entries with no MatchedFragmentIons
        /// </summary>
        /// <returns>Array of LibraryPrecursorInput.</returns>
        public static LibraryPrecursorInput[] ConvertLibrarySpectra(
            List<LibrarySpectrum> librarySpectra)
        {
            if (librarySpectra == null || librarySpectra.Count == 0)
                return Array.Empty<LibraryPrecursorInput>();

            var inputs = new List<LibraryPrecursorInput>(librarySpectra.Count);

            for (int i = 0; i < librarySpectra.Count; i++)
            {
                var lib = librarySpectra[i];

                // Skip entries with no fragment ions
                if (lib.MatchedFragmentIons == null || lib.MatchedFragmentIons.Count == 0)
                    continue;

                int fragCount = lib.MatchedFragmentIons.Count;
                var fragMzs = new float[fragCount];
                var fragIntensities = new float[fragCount];

                for (int f = 0; f < fragCount; f++)
                {
                    fragMzs[f] = (float)lib.MatchedFragmentIons[f].Mz;
                    fragIntensities[f] = (float)lib.MatchedFragmentIons[f].Intensity;
                }

                inputs.Add(new LibraryPrecursorInput(
                    sequence: lib.Sequence,
                    precursorMz: lib.PrecursorMz,           // double — precision preserved
                    chargeState: lib.ChargeState,
                    retentionTime: lib.RetentionTime,        // double? — precision preserved
                    isDecoy: lib.IsDecoy,
                    fragmentMzs: fragMzs,                    // float[] — converted from double
                    fragmentIntensities: fragIntensities));   // float[] — converted from double
            }

            return inputs.ToArray();
        }

        /// <summary>
        /// Disposes the DIA scan index if still held.
        /// Normally cleaned up at the end of RunSpecific(), but this provides
        /// safety in case of early termination.
        /// </summary>
        public void Dispose()
        {
            if (!_disposed)
            {
                _index?.Dispose();
                _index = null;
                _disposed = true;
            }
        }
    }
}