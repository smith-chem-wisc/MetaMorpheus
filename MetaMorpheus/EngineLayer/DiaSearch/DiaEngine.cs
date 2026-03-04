// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using MassSpectrometry;
using MassSpectrometry.Dia;
using Omics.SpectrumMatch;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.DiaSearch
{
    /// <summary>
    /// MetaMorpheus engine for DIA library-based spectral matching with iRT calibration.
    /// 
    /// Orchestrates the full DIA search pipeline:
    ///   Phase 0: Build SoA scan index from MsDataFile, convert LibrarySpectrum → LibraryPrecursorInput
    ///   Phase 1: Broad search with provisional iRT mapping (wide windows)
    ///   Phase 1b: Extract anchor matches → fit iRT calibration model (RANSAC + OLS)
    ///   Phase 2+: Refined search with calibrated RT windows (k·σ_iRT) and RT scoring
    ///   Final: Assemble DiaSearchResult list with combined spectral + RT scores
    /// 
    /// The engine delegates all heavy computation to mzLib components:
    ///   - DiaScanIndexBuilder: SoA data ingestion
    ///   - DiaExtractionOrchestrator: parallel window-level fragment extraction
    ///   - DiaLibraryQueryGenerator: query generation, anchor extraction, result assembly
    ///   - RtCalibrationFitter: robust linear regression
    ///   - NormalizedDotProductScorer / SpectralAngleScorer: spectral scoring
    /// 
    /// This class is the bridge between MetaMorpheus's Omics types (LibrarySpectrum,
    /// MatchedFragmentIon) and mzLib's DIA engine types (LibraryPrecursorInput, FragmentQuery).
    /// </summary>
    public class DiaEngine : MetaMorpheusEngine, IDisposable
    {
        private readonly MsDataScan[] _scans;
        private readonly List<LibrarySpectrum> _librarySpectra;
        private readonly DiaSearchParameters _diaParameters;
        private readonly string _fileName;

        /// <summary>
        /// Creates a DIA search engine.
        /// </summary>
        /// <param name="dataFile">The DIA raw data file (mzML, Thermo raw, etc.).</param>
        /// <param name="librarySpectra">
        /// Library spectra from Koina/Prosit predictions or .msp files.
        /// Each spectrum's RetentionTime is treated as iRT when UseIrtCalibration is enabled.
        /// </param>
        /// <param name="diaParameters">DIA-specific search parameters.</param>
        /// <param name="commonParameters">MetaMorpheus common parameters.</param>
        /// <param name="fileSpecificParameters">File-specific parameter overrides.</param>
        /// <param name="nestedIds">Engine nesting identifiers for progress reporting.</param>
        public DiaEngine(
            MsDataScan[] scans,
            List<LibrarySpectrum> librarySpectra,
            DiaSearchParameters diaParameters,
            CommonParameters commonParameters,
            List<(string FileName, CommonParameters Parameters)> fileSpecificParameters,
            List<string> nestedIds,
            string fileName = "unknown")
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            _scans = scans ?? throw new ArgumentNullException(nameof(scans));
            _librarySpectra = librarySpectra ?? throw new ArgumentNullException(nameof(librarySpectra));
            _diaParameters = diaParameters ?? throw new ArgumentNullException(nameof(diaParameters));
            _fileName = fileName;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            // ── Phase 0: Build scan index and convert library ──────────────────
            Status("Building DIA scan index...");
            var scanIndex = DiaScanIndexBuilder.Build(_scans);

            if (scanIndex.ScanCount == 0)
            {
                Warn("No MS2 scans found in DIA file. Returning empty results.");
                return new DiaEngineResults(this, new List<DiaSearchResult>(), null, 0, 0);
            }

            Status($"DIA index: {scanIndex.ScanCount} scans, {scanIndex.TotalPeakCount:N0} peaks, {scanIndex.WindowCount} windows");

            // Convert LibrarySpectrum → LibraryPrecursorInput
            Status("Converting library spectra...");
            var precursors = ConvertLibrarySpectra(_librarySpectra, _diaParameters.UseIrtCalibration);

            // Count targets and decoys from the library (before any filtering by the search)
            int targetLibraryCount = precursors.Count(p => !p.IsDecoy);
            int decoyLibraryCount = precursors.Count(p => p.IsDecoy);

            if (precursors.Count == 0)
            {
                Warn("No valid library precursors after conversion. Returning empty results.");
                return new DiaEngineResults(this, new List<DiaSearchResult>(), null, 0, 0);
            }

            Status($"Library: {precursors.Count} precursors ({targetLibraryCount} targets, {decoyLibraryCount} decoys)");

            // Create scorers
            var dotProductScorer = new NormalizedDotProductScorer();
            var spectralAngleScorer = new SpectralAngleScorer();

            // Create orchestrator
            var orchestrator = new DiaExtractionOrchestrator(
                scanIndex,
                _diaParameters.PreferGpu
                    ? FragmentExtractorFactory.CreateFactory()
                    : (idx => new CpuFragmentExtractor(idx)));

            int maxThreads = _diaParameters.EffectiveMaxThreads;
            RtCalibrationModel calibration = null;
            List<DiaSearchResult> results = null;

            try
            {
                if (_diaParameters.UseIrtCalibration)
                {
                    // ── iRT Calibration Pipeline ──────────────────────────────────
                    results = RunCalibratedSearch(
                        precursors, scanIndex, orchestrator, maxThreads,
                        dotProductScorer, spectralAngleScorer,
                        out calibration);
                }
                else
                {
                    // ── Fixed RT Window Pipeline ──────────────────────────────────
                    results = RunFixedRtSearch(
                        precursors, scanIndex, orchestrator, maxThreads,
                        dotProductScorer, spectralAngleScorer);
                }
            }
            finally
            {
                (orchestrator as IDisposable)?.Dispose();
                (scanIndex as IDisposable)?.Dispose();
            }

            Status($"DIA search complete: {results.Count} results " +
                   $"({results.Count(r => !r.IsDecoy)} targets, {results.Count(r => r.IsDecoy)} decoys)");

            return new DiaEngineResults(this, results, calibration, targetLibraryCount, decoyLibraryCount);
        }

        // ── Two-Phase Calibrated Search ─────────────────────────────────────────

        private List<DiaSearchResult> RunCalibratedSearch(
            List<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            DiaExtractionOrchestrator orchestrator,
            int maxThreads,
            IScorer dotProductScorer,
            IScorer spectralAngleScorer,
            out RtCalibrationModel calibration)
        {
            // Build iRT library index for candidate selection
            var irtIndex = new IrtLibraryIndex(precursors);

            Status($"Run RT range: {scanIndex.GetGlobalRtMin():F2} to {scanIndex.GetGlobalRtMax():F2} minutes");

            if (irtIndex.Count == 0)
            {
                Warn("No precursors have iRT values. Falling back to fixed RT window search.");
                calibration = null;
                return RunFixedRtSearch(precursors, scanIndex, orchestrator, maxThreads,
                    dotProductScorer, spectralAngleScorer);
            }

            // Phase 1: Provisional calibration (broad windows)
            Status("Phase 1: Broad search with provisional iRT mapping...");
            var provisionalCalibration = RtCalibrationModel.CreateProvisional(
                runRtMin: scanIndex.GetGlobalRtMin(),
                runRtMax: scanIndex.GetGlobalRtMax(),
                libraryIrtMin: irtIndex.MinIrt,
                libraryIrtMax: irtIndex.MaxIrt,
                initialWindowIrt: _diaParameters.InitialIrtWindow);

            Status($"Provisional calibration: {provisionalCalibration}");


            // Generate queries with provisional calibration
            var genResult = DiaLibraryQueryGenerator.GenerateCalibrated(
                precursors, scanIndex, _diaParameters, provisionalCalibration);

            Status($"Phase 1 queries: {genResult.Queries.Length} " +
                   $"({genResult.PrecursorGroups.Length} precursors, " +
                   $"skipped: {genResult.SkippedNoWindow} no-window, {genResult.SkippedNoFragments} no-fragments)");


            Warn($"DEBUG Run RT range: {scanIndex.GetGlobalRtMin():F2} to {scanIndex.GetGlobalRtMax():F2} minutes");

            // Extract
            var extractionResult = orchestrator.ExtractAll(genResult.Queries, maxThreads);

            // Assemble Phase 1 results (no RT scoring yet — provisional calibration is unreliable)
            var phase1Results = DiaLibraryQueryGenerator.AssembleResults(
                precursors, genResult, extractionResult, _diaParameters,
                dotProductScorer, spectralAngleScorer, calibration: null);

            Status($"Phase 1 results: {phase1Results.Count}");

            //debug
            // Diagnostic: score distribution of Phase 1 results
            if (phase1Results.Count > 0)
            {
                var scores = phase1Results.Select(r => r.DotProductScore).Where(s => !float.IsNaN(s)).OrderByDescending(s => s).ToList();
                Status($"Phase 1 score distribution: max={scores.FirstOrDefault():F3}, " +
                       $"median={scores[scores.Count / 2]:F3}, " +
                       $"≥0.5: {scores.Count(s => s >= 0.5)}, " +
                       $"≥0.3: {scores.Count(s => s >= 0.3)}, " +
                       $"≥0.1: {scores.Count(s => s >= 0.1)}");

                // Check how many targets vs decoys have iRT values
                var targetsWithIrt = phase1Results.Count(r => !r.IsDecoy && r.LibraryRetentionTime.HasValue);
                Status($"Phase 1 targets with RT: {targetsWithIrt}/{phase1Results.Count(r => !r.IsDecoy)}");
            }


            // Extract anchors from Phase 1
            DiaLibraryQueryGenerator.ExtractAnchors(
                phase1Results,
                _diaParameters.CalibrationAnchorMinScore,
                out double[] anchorIrts,
                out double[] anchorRtMinutes);

            Status($"Calibration anchors: {anchorIrts.Length}");

            if (anchorIrts.Length < RtCalibrationModel.MinReliableAnchors)
            {
                Warn($"Only {anchorIrts.Length} anchors found (minimum: {RtCalibrationModel.MinReliableAnchors}). " +
                     "Falling back to fixed RT window search.");
                calibration = null;
                return RunFixedRtSearch(precursors, scanIndex, orchestrator, maxThreads,
                    dotProductScorer, spectralAngleScorer);
            }

            // Fit calibration
            calibration = RtCalibrationFitter.Fit((IList<double>)anchorIrts, (IList<double>)anchorRtMinutes);

            if (calibration == null || !calibration.IsReliable)
            {
                Warn($"Calibration fitting failed or unreliable " +
                     $"(R²={(calibration?.RSquared ?? 0):F4}, anchors={calibration?.AnchorCount ?? 0}). " +
                     "Falling back to Phase 1 results.");
                calibration = null;
                return phase1Results;
            }

            Status($"Calibration fit: {calibration}");

            // Phase 2+: Iterative refinement
            List<DiaSearchResult> currentResults = phase1Results;
            double prevSlope = calibration.Slope;

            for (int iteration = 1; iteration <= _diaParameters.MaxCalibrationIterations; iteration++)
            {
                Status($"Calibration iteration {iteration}: re-generating queries with calibrated windows...");

                // Re-generate queries with calibrated windows
                genResult = DiaLibraryQueryGenerator.GenerateCalibrated(
                    precursors, scanIndex, _diaParameters, calibration);

                Status($"Iteration {iteration} queries: {genResult.Queries.Length}");

                // Re-extract
                extractionResult = orchestrator.ExtractAll(genResult.Queries, maxThreads);

                // Assemble with RT scoring
                currentResults = DiaLibraryQueryGenerator.AssembleResults(
                    precursors, genResult, extractionResult, _diaParameters,
                    dotProductScorer, spectralAngleScorer, calibration);

                Status($"Iteration {iteration} results: {currentResults.Count}");

                // Re-extract anchors and refit
                DiaLibraryQueryGenerator.ExtractAnchors(
                    currentResults,
                    _diaParameters.CalibrationAnchorMinScore,
                    out anchorIrts,
                    out anchorRtMinutes);

                Status($"Iteration {iteration} anchors: {anchorIrts.Length}");

                if (anchorIrts.Length < RtCalibrationModel.MinReliableAnchors)
                    break;

                var newCalibration = RtCalibrationFitter.Fit((IList<double>)anchorIrts, (IList<double>)anchorRtMinutes);

                if (newCalibration == null || !newCalibration.IsReliable)
                    break;

                // Check convergence
                double slopeDelta = Math.Abs(newCalibration.Slope - prevSlope);
                prevSlope = newCalibration.Slope;
                calibration = newCalibration;

                Status($"Iteration {iteration} calibration: {calibration} (Δslope={slopeDelta:E3})");

                if (slopeDelta < _diaParameters.CalibrationConvergenceEpsilon)
                {
                    Status($"Calibration converged after {iteration} iterations.");
                    break;
                }
            }

            // Final pass with converged calibration (if the last iteration didn't already produce it)
            if (currentResults == phase1Results)
            {
                // We never successfully ran a calibrated pass — use phase1 results
                return currentResults;
            }

            return currentResults;
        }

        // ── Fixed RT Window Search (fallback) ───────────────────────────────────

        private List<DiaSearchResult> RunFixedRtSearch(
            List<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            DiaExtractionOrchestrator orchestrator,
            int maxThreads,
            IScorer dotProductScorer,
            IScorer spectralAngleScorer)
        {
            Status("Running fixed RT window search (no iRT calibration)...");

            var genResult = DiaLibraryQueryGenerator.Generate(
                precursors, scanIndex, _diaParameters);

            Status($"Queries: {genResult.Queries.Length} " +
                   $"({genResult.PrecursorGroups.Length} precursors, " +
                   $"skipped: {genResult.SkippedNoWindow} no-window, {genResult.SkippedNoFragments} no-fragments)");

            var extractionResult = orchestrator.ExtractAll(genResult.Queries, maxThreads);

            var results = DiaLibraryQueryGenerator.AssembleResults(
                precursors, genResult, extractionResult, _diaParameters,
                dotProductScorer, spectralAngleScorer);

            return results;
        }

        // ── Library Conversion ──────────────────────────────────────────────────

        /// <summary>
        /// Converts MetaMorpheus LibrarySpectrum objects to mzLib's LibraryPrecursorInput structs.
        /// 
        /// This is the bridge between Omics types and the MassSpectrometry.Dia engine.
        /// When iRT calibration is enabled, LibrarySpectrum.RetentionTime is treated as iRT
        /// (since Koina/Prosit predictions output iRT values in the RetentionTime field).
        /// </summary>
        internal static List<LibraryPrecursorInput> ConvertLibrarySpectra(
            IList<LibrarySpectrum> spectra,
            bool treatRetentionTimeAsIrt = true)
        {
            if (spectra == null)
                return new List<LibraryPrecursorInput>();

            var precursors = new List<LibraryPrecursorInput>(spectra.Count);

            for (int i = 0; i < spectra.Count; i++)
            {
                var lib = spectra[i];

                if (lib.MatchedFragmentIons == null || lib.MatchedFragmentIons.Count == 0)
                    continue;

                // Extract fragment m/z and intensity arrays
                int fragCount = lib.MatchedFragmentIons.Count;
                var mzs = new float[fragCount];
                var intensities = new float[fragCount];

                for (int f = 0; f < fragCount; f++)
                {
                    var ion = lib.MatchedFragmentIons[f];
                    mzs[f] = (float)ion.Mz;
                    intensities[f] = (float)ion.Intensity;
                }

                // When iRT calibration is active, the library's RetentionTime IS the iRT
                // (Koina/Prosit outputs iRT in the RetentionTime field)
                double? irtValue = treatRetentionTimeAsIrt ? lib.RetentionTime : null;

                precursors.Add(new LibraryPrecursorInput(
                    sequence: lib.Sequence,
                    precursorMz: lib.PrecursorMz,
                    chargeState: lib.ChargeState,
                    retentionTime: lib.RetentionTime,
                    isDecoy: lib.IsDecoy,
                    fragmentMzs: mzs,
                    fragmentIntensities: intensities,
                    irtValue: irtValue));
            }

            return precursors;
        }
        public void Dispose()
        {
            // DiaEngine does not own any unmanaged resources.
            // The DiaScanIndex and DiaExtractionOrchestrator are created and disposed
            // within RunSpecific(), so nothing to clean up here.
            // This interface is required because callers use `using` statements.
        }
    }
}