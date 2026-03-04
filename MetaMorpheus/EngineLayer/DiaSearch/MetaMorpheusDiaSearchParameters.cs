// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using MassSpectrometry.Dia;

namespace EngineLayer.DiaSearch
{
    /// <summary>
    /// MetaMorpheus-side wrapper for DIA search parameters.
    /// Exposes settings for TOML serialization and provides a factory
    /// method to create the mzLib DiaSearchParameters used by the engine.
    /// 
    /// This class exists because DiaSearchParameters lives in mzLib (MassSpectrometry.Dia)
    /// which has no knowledge of MetaMorpheus's TOML config system. This wrapper adds
    /// the serialization layer and any MetaMorpheus-specific defaults.
    /// </summary>
    public class MetaMorpheusDiaSearchParameters
    {
        // ── Fragment Extraction ──────────────────────────────────────────────

        /// <summary>Fragment m/z tolerance in ppm. Default: 20.</summary>
        public float PpmTolerance { get; set; } = 20f;

        /// <summary>
        /// Fallback RT window half-width in minutes (used when iRT calibration is disabled or fails).
        /// Default: 5.0.
        /// </summary>
        public float RtToleranceMinutes { get; set; } = 5.0f;

        /// <summary>Minimum detected fragments for a result to be reported. Default: 3.</summary>
        public int MinFragmentsRequired { get; set; } = 3;

        /// <summary>Minimum score threshold (pre-FDR). Default: 0.0.</summary>
        public float MinScoreThreshold { get; set; } = 0.0f;

        // ── Parallelism ─────────────────────────────────────────────────────

        /// <summary>Max threads for extraction. -1 = all cores. Default: -1.</summary>
        public int MaxThreads { get; set; } = -1;

        /// <summary>Whether to attempt GPU acceleration. Default: false.</summary>
        public bool PreferGpu { get; set; } = false;

        // ── iRT Calibration ─────────────────────────────────────────────────

        /// <summary>
        /// Enable iRT-based retention time calibration and scoring.
        /// When true, LibrarySpectrum.RetentionTime is treated as iRT.
        /// Default: true.
        /// </summary>
        public bool UseIrtCalibration { get; set; } = true;

        /// <summary>
        /// Broad window half-width in iRT units for Phase 1 (before calibration).
        /// Default: 20.
        /// </summary>
        public double InitialIrtWindow { get; set; } = 40.0;

        /// <summary>
        /// Number of σ for calibrated RT window: window = k · σ_iRT.
        /// Default: 3.0.
        /// </summary>
        public double CalibratedWindowSigmaMultiplier { get; set; } = 3.0;

        /// <summary>
        /// Minimum spectral score to qualify as a calibration anchor.
        /// Default: 0.5.
        /// </summary>
        public float CalibrationAnchorMinScore { get; set; } = 0.3f;

        /// <summary>
        /// Weight (λ) for RT score in combined score: finalScore = spectral + λ · rtScore.
        /// 0 disables RT scoring. Default: 0.5.
        /// </summary>
        public double RtScoreLambda { get; set; } = 0.5;

        /// <summary>
        /// Maximum calibration refinement iterations. Default: 3.
        /// </summary>
        public int MaxCalibrationIterations { get; set; } = 3;

        /// <summary>
        /// Slope convergence threshold for stopping iteration. Default: 0.001.
        /// </summary>
        public double CalibrationConvergenceEpsilon { get; set; } = 0.001;

        // ── Output Control ──────────────────────────────────────────────────
        // These are MetaMorpheus-specific (not passed to mzLib). They control
        // what gets written to disk by DiaSearchTask after the engine runs.
        // Follows the same pattern as SearchParameters.WriteDecoys, etc.

        /// <summary>
        /// Write diagnostic output files (calibration plots, per-window statistics,
        /// anchor lists, score distributions). Useful for debugging and method development.
        /// Default: false.
        /// </summary>
        public bool WriteDiagnostics { get; set; } = true;

        /// <summary>
        /// Include decoy results in the output TSV files.
        /// When false, only target identifications are written.
        /// Decoys are always used internally for FDR estimation regardless of this setting.
        /// Default: true (matching SearchParameters.WriteDecoys convention).
        /// </summary>
        public bool WriteDecoyResults { get; set; } = false;

        /// <summary>
        /// Include contaminant results in the output TSV files.
        /// Default: true (matching SearchParameters.WriteContaminants convention).
        /// </summary>
        public bool WriteContaminantResults { get; set; } = true;

        /// <summary>
        /// Write per-file result files in addition to the combined results.
        /// Default: true (matching SearchParameters.WriteIndividualFiles convention).
        /// </summary>
        public bool WriteIndividualFiles { get; set; } = true;

        /// <summary>
        /// Write results that exceed the q-value threshold (high q-value results).
        /// Useful for inspection but typically filtered out in final analysis.
        /// Default: true (matching SearchParameters.WriteHighQValuePsms convention).
        /// </summary>
        public bool WriteHighQValueResults { get; set; } = true;
        /// <summary>
        /// Path to a spectral library file (.msp) to use for DIA searching.
        /// When null/empty, library spectra are generated from Koina/Prosit predictions.
        /// </summary>
        public string SpectralLibraryPath { get; set; } = "";
        /// <summary>
        /// Creates the mzLib DiaSearchParameters from this configuration.
        /// Only passes engine-relevant parameters; output control flags stay in MetaMorpheus.
        /// </summary>
        public DiaSearchParameters ToMzLibParameters()
        {
            return new DiaSearchParameters
            {
                PpmTolerance = PpmTolerance,
                RtToleranceMinutes = RtToleranceMinutes,
                MinFragmentsRequired = MinFragmentsRequired,
                MinScoreThreshold = MinScoreThreshold,
                MaxThreads = MaxThreads,
                PreferGpu = PreferGpu,
                UseIrtCalibration = UseIrtCalibration,
                InitialIrtWindow = InitialIrtWindow,
                CalibratedWindowSigmaMultiplier = CalibratedWindowSigmaMultiplier,
                CalibrationAnchorMinScore = CalibrationAnchorMinScore,
                RtScoreLambda = RtScoreLambda,
                MaxCalibrationIterations = MaxCalibrationIterations,
                CalibrationConvergenceEpsilon = CalibrationConvergenceEpsilon
            };
        }
    }
}