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
    /// MetaMorpheus engine for DIA library-based spectral matching with iterative RT calibration.
    ///
    /// Delegates all heavy computation to mzLib via two entry points:
    ///
    ///   Calibrated path (default, useCalibration=true):
    ///     DiaCalibrationPipeline.RunWithAutomaticCalibration — runs the full bootstrap →
    ///     iterative calibration → final adaptive extraction → temporal scoring pipeline
    ///     in a single call. This is the path that produces 29,157 IDs at 1% FDR.
    ///
    ///   Fixed RT window path (useCalibration=false):
    ///     DiaLibraryQueryGenerator.Generate + ExtractAll + AssembleResults — no calibration.
    ///
    /// The scan index is deliberately NOT disposed here — it is passed out through
    /// DiaEngineResults.ScanIndex so PostDiaSearchAnalysisTask can supply it to
    /// DiaFeatureExtractor for MS1 feature computation (isotope scores, MS1-MS2
    /// correlation, etc.). DiaSearchTask disposes it after FDR completes.
    /// </summary>
    public class DiaEngine : MetaMorpheusEngine, IDisposable
    {
        private readonly MsDataScan[] _scans;
        private readonly List<LibrarySpectrum> _librarySpectra;
        private readonly DiaSearchParameters _diaParameters;
        private readonly bool _useCalibration;
        private readonly string _fileName;

        /// <param name="scans">All scans from the DIA raw file.</param>
        /// <param name="librarySpectra">Target + decoy library spectra (IsDecoy must be set).</param>
        /// <param name="diaParameters">mzLib DIA parameters.</param>
        /// <param name="useCalibration">
        /// True (default) = use DiaCalibrationPipeline for iterative RT calibration.
        /// False = fixed RT window fallback (faster, fewer IDs).
        /// </param>
        /// <param name="commonParameters">MetaMorpheus common parameters.</param>
        /// <param name="fileSpecificParameters">File-specific parameter overrides.</param>
        /// <param name="nestedIds">Progress reporting identifiers.</param>
        /// <param name="fileName">Raw file name for log messages.</param>
        public DiaEngine(
            MsDataScan[] scans,
            List<LibrarySpectrum> librarySpectra,
            DiaSearchParameters diaParameters,
            bool useCalibration,
            CommonParameters commonParameters,
            List<(string FileName, CommonParameters Parameters)> fileSpecificParameters,
            List<string> nestedIds,
            string fileName = "unknown")
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            _scans = scans ?? throw new ArgumentNullException(nameof(scans));
            _librarySpectra = librarySpectra ?? throw new ArgumentNullException(nameof(librarySpectra));
            _diaParameters = diaParameters ?? throw new ArgumentNullException(nameof(diaParameters));
            _useCalibration = useCalibration;
            _fileName = fileName;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            // ── Phase 0: Build scan index ────────────────────────────────────────
            Status("Building DIA scan index...");
            var scanIndex = DiaScanIndexBuilder.Build(_scans);

            if (scanIndex.ScanCount == 0)
            {
                Warn("No MS2 scans found in DIA file. Returning empty results.");
                return new DiaEngineResults(this, new List<DiaSearchResult>(), null, null, 0, 0);
            }

            Status($"DIA index: {scanIndex.ScanCount} scans, {scanIndex.TotalPeakCount:N0} peaks, " +
                   $"{scanIndex.WindowCount} windows");

            // ── Phase 0b: Convert library types ──────────────────────────────────
            Status("Converting library spectra...");
            var precursors = ConvertLibrarySpectra(_librarySpectra);
            int targetLibraryCount = precursors.Count(p => !p.IsDecoy);
            int decoyLibraryCount = precursors.Count(p => p.IsDecoy);

            if (precursors.Count == 0)
            {
                Warn("No valid library precursors after conversion. Returning empty results.");
                return new DiaEngineResults(this, new List<DiaSearchResult>(), null, null, 0, 0);
            }

            Status($"Library: {precursors.Count} precursors " +
                   $"({targetLibraryCount} targets, {decoyLibraryCount} decoys)");

            // ── Phase 1+: Search and calibration ─────────────────────────────────
            RtCalibrationModel calibration = null;
            List<DiaSearchResult> results;

            var orchestrator = new DiaExtractionOrchestrator(
                scanIndex,
                _diaParameters.PreferGpu
                    ? FragmentExtractorFactory.CreateFactory()
                    : (idx => new CpuFragmentExtractor(idx)));

            try
            {
                if (_useCalibration)
                {
                    System.Diagnostics.Debug.WriteLine("Running calibrated DIA search (iterative RT calibration)...");

                    var pipelineResult = DiaCalibrationPipeline.RunWithAutomaticCalibration(
                        precursors, scanIndex, _diaParameters, orchestrator);

                    results = pipelineResult.Results;
                    calibration = pipelineResult.Calibration;

                    // ── DEBUG: confirm calibration ran and what it produced ─────────────
                    if (pipelineResult.CalibrationLog != null)
                    {
                        System.Diagnostics.Debug.WriteLine($"[DIA DEBUG] Calibration log entries: {pipelineResult.CalibrationLog.Count}");
                        foreach (var entry in pipelineResult.CalibrationLog)
                            System.Diagnostics.Debug.WriteLine($"[DIA DEBUG] Iter: sigma={entry.SigmaMinutes:F4} R2={entry.RSquared:F4} window={entry.WindowHalfWidthMinutes:F4}");
                    }
                    else
                    {
                        System.Diagnostics.Debug.WriteLine("[DIA DEBUG] CalibrationLog is NULL");
                    }
                    System.Diagnostics.Debug.WriteLine($"[DIA DEBUG] Calibration model null={pipelineResult.Calibration == null}");
                    System.Diagnostics.Debug.WriteLine($"[DIA DEBUG] Results count={pipelineResult.Results?.Count ?? -1}");
                }
                else
                {
                    Status("Running fixed RT window DIA search (no calibration)...");
                    results = RunFixedRtSearch(precursors, scanIndex, orchestrator);
                }
            }
            finally
            {
                // Dispose orchestrator but NOT scanIndex — it flows out via DiaEngineResults
                // and is needed for MS1 feature extraction in PostDiaSearchAnalysisTask.
                (orchestrator as IDisposable)?.Dispose();
            }

            Status($"DIA search complete: {results.Count} results " +
                   $"({results.Count(r => !r.IsDecoy)} targets, {results.Count(r => r.IsDecoy)} decoys)");

            return new DiaEngineResults(
                this, results, calibration, scanIndex,
                targetLibraryCount, decoyLibraryCount);
        }

        // ── Fixed RT Window Search ──────────────────────────────────────────────

        private List<DiaSearchResult> RunFixedRtSearch(
            List<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            DiaExtractionOrchestrator orchestrator)
        {
            var genResult = DiaLibraryQueryGenerator.Generate(
                precursors, scanIndex, _diaParameters);

            Status($"Fixed RT queries: {genResult.Queries.Length} " +
                   $"({genResult.PrecursorGroups.Length} precursors, " +
                   $"skipped: {genResult.SkippedNoWindow} no-window, " +
                   $"{genResult.SkippedNoFragments} no-fragments)");

            // ExtractAll returns ExtractionResult — pass .Results (FragmentResult[]) to AssembleResults
            var extractionResult = orchestrator.ExtractAll(
                genResult.Queries, _diaParameters.EffectiveMaxThreads);

            return DiaLibraryQueryGenerator.AssembleResults(
                precursors,
                genResult,
                extractionResult.Results,
                _diaParameters,
                dotProductScorer: new NormalizedDotProductScorer(),
                spectralAngleScorer: new SpectralAngleScorer());
        }

        // ── Library Conversion ──────────────────────────────────────────────────

        /// <summary>
        /// Converts MetaMorpheus LibrarySpectrum objects to mzLib LibraryPrecursorInput.
        ///
        /// LibrarySpectrum.RetentionTime is passed through as-is. When the library comes
        /// from Koina/Prosit predictions, this field contains the predicted iRT value and
        /// the DiaCalibrationPipeline uses it for RT calibration.
        ///
        /// Note: LibraryPrecursorInput has no iRT-specific field — the calibration pipeline
        /// infers the RT-to-iRT relationship from the distribution of library RetentionTime
        /// values vs observed apex RTs.
        /// </summary>
        internal static List<LibraryPrecursorInput> ConvertLibrarySpectra(
            IList<LibrarySpectrum> spectra)
        {
            if (spectra == null) return new List<LibraryPrecursorInput>();

            var precursors = new List<LibraryPrecursorInput>(spectra.Count);

            for (int i = 0; i < spectra.Count; i++)
            {
                var lib = spectra[i];

                if (lib.MatchedFragmentIons == null || lib.MatchedFragmentIons.Count == 0)
                    continue;

                int fragCount = lib.MatchedFragmentIons.Count;
                var mzs = new float[fragCount];
                var intensities = new float[fragCount];

                for (int f = 0; f < fragCount; f++)
                {
                    mzs[f] = (float)lib.MatchedFragmentIons[f].Mz;
                    intensities[f] = (float)lib.MatchedFragmentIons[f].Intensity;
                }

                precursors.Add(new LibraryPrecursorInput(
                    sequence: lib.Sequence,
                    precursorMz: lib.PrecursorMz,
                    chargeState: lib.ChargeState,
                    retentionTime: lib.RetentionTime,
                    isDecoy: lib.IsDecoy,
                    fragmentMzs: mzs,
                    fragmentIntensities: intensities));
            }

            return precursors;
        }

        public void Dispose() { }
    }
}