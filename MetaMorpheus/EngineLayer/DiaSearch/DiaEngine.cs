// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using MassSpectrometry;
using MassSpectrometry.Dia;
using Omics.Fragmentation;
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
    /// Scorers).
    /// 
    /// Pipeline (executed in RunSpecific):
    ///   1. Convert LibrarySpectrum list → LibraryPrecursorInput[]
    ///   2. Build DiaScanIndex from MsDataScan[]
    ///   3. Generate FragmentQuery[] via DiaLibraryQueryGenerator
    ///   4. Extract via DiaExtractionOrchestrator (parallel, CPU or GPU)
    ///   5. Score and assemble DiaSearchResult list
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
            Status("Converting spectral library...");
            var precursorInputs = ConvertLibrarySpectra(_librarySpectra);

            int targetCount = 0;
            int decoyCount = 0;
            for (int i = 0; i < precursorInputs.Length; i++)
            {
                if (precursorInputs[i].IsDecoy) decoyCount++;
                else targetCount++;
            }

            if (precursorInputs.Length == 0)
            {
                return new DiaEngineResults(this,
                    new List<DiaSearchResult>(),
                    targetLibraryCount: 0,
                    decoyLibraryCount: 0);
            }

            // ── Step 1: Build SoA index ─────────────────────────────────────
            Status("Building DIA scan index...");
            _index = DiaScanIndexBuilder.Build(_scans);

            Status($"Indexed {_index.ScanCount} MS2 scans across {_index.WindowCount} windows " +
                   $"({_index.TotalPeakCount:N0} peaks)");

            // ── Step 2: Generate queries ────────────────────────────────────
            Status("Generating fragment queries...");
            var generationResult = DiaLibraryQueryGenerator.Generate(
                precursorInputs, _index, _diaParameters);

            Status($"Generated {generationResult.Queries.Length:N0} queries " +
                   $"({generationResult.SkippedNoWindow} precursors outside windows)");

            if (generationResult.Queries.Length == 0)
            {
                _index.Dispose();
                return new DiaEngineResults(this,
                    new List<DiaSearchResult>(),
                    targetLibraryCount: targetCount,
                    decoyLibraryCount: decoyCount);
            }

            // ── Step 3: Extract fragment XICs ───────────────────────────────
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

            Status($"DIA search complete: {results.Count:N0} precursor matches " +
                   $"({results.Count(r => !r.IsDecoy)} targets, {results.Count(r => r.IsDecoy)} decoys)");

            _index.Dispose();

            return new DiaEngineResults(this, results,
                targetLibraryCount: targetCount,
                decoyLibraryCount: decoyCount);
        }

        /// <summary>
        /// Converts Omics LibrarySpectrum objects into mzLib-native LibraryPrecursorInput
        /// structs that have no Omics dependency.
        /// </summary>
        /// <returns>Array of LibraryPrecursorInput (skips entries with no fragments).</returns>
        public static LibraryPrecursorInput[] ConvertLibrarySpectra(
            List<LibrarySpectrum> librarySpectra)
        {
            if (librarySpectra == null || librarySpectra.Count == 0)
                return Array.Empty<LibraryPrecursorInput>();

            var inputs = new List<LibraryPrecursorInput>(librarySpectra.Count);

            for (int i = 0; i < librarySpectra.Count; i++)
            {
                var lib = librarySpectra[i];
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
                    precursorMz: lib.PrecursorMz,
                    chargeState: lib.ChargeState,
                    retentionTime: lib.RetentionTime,
                    isDecoy: lib.IsDecoy,
                    fragmentMzs: fragMzs,
                    fragmentIntensities: fragIntensities));
            }

            return inputs.ToArray();
        }

        public void Dispose()
        {
            if (!_disposed)
            {
                _index?.Dispose();
                _disposed = true;
            }
        }
    }

    /// <summary>
    /// Results container for DiaEngine, extending MetaMorpheusEngineResults
    /// for integration with the standard MetaMorpheus reporting pipeline.
    /// </summary>
    public sealed class DiaEngineResults : MetaMorpheusEngineResults
    {
        /// <summary>DIA precursor match results (targets + decoys that passed filters).</summary>
        public List<DiaSearchResult> DiaResults { get; }

        /// <summary>Number of target (non-decoy) library entries provided as input.</summary>
        public int TargetLibraryCount { get; }

        /// <summary>Number of decoy library entries provided as input.</summary>
        public int DecoyLibraryCount { get; }

        public DiaEngineResults(
            DiaEngine engine,
            List<DiaSearchResult> diaResults,
            int targetLibraryCount,
            int decoyLibraryCount)
            : base(engine)
        {
            DiaResults = diaResults ?? new List<DiaSearchResult>();
            TargetLibraryCount = targetLibraryCount;
            DecoyLibraryCount = decoyLibraryCount;
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine($"Library: {TargetLibraryCount} targets, {DecoyLibraryCount} decoys");
            sb.AppendLine($"DIA matches: {DiaResults.Count}");

            if (DiaResults.Count > 0)
            {
                int targetHits = DiaResults.Count(r => !r.IsDecoy);
                int decoyHits = DiaResults.Count(r => r.IsDecoy);
                sb.AppendLine($"  Target hits: {targetHits}, Decoy hits: {decoyHits}");

                var scored = DiaResults.Where(r => !float.IsNaN(r.DotProductScore)).ToList();
                if (scored.Count > 0)
                {
                    sb.AppendLine($"  Avg dot product: {scored.Average(r => r.DotProductScore):F4}");
                    sb.AppendLine($"  Best dot product: {scored.Max(r => r.DotProductScore):F4}");
                }
            }

            return sb.ToString();
        }
    }
}
