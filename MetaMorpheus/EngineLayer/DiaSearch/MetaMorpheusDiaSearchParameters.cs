// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using MassSpectrometry.Dia;

namespace EngineLayer.DiaSearch
{
    /// <summary>
    /// MetaMorpheus-side configuration for DIA search.
    ///
    /// Exposes settings for TOML serialization and GUI binding, and provides
    /// ToMzLibParameters() to create the mzLib DiaSearchParameters used by DiaEngine.
    ///
    /// Design note: DiaSearchParameters (mzLib) intentionally has no RT calibration
    /// fields — calibration is fully automatic and controlled by IterativeRtCalibrator
    /// inside DiaCalibrationPipeline. The only RT parameter exposed here is
    /// RtToleranceMinutes, which sets the initial broad-window half-width for the
    /// bootstrap pass (iteration 0). After iteration 0, adaptive windows are used.
    /// </summary>
    public class MetaMorpheusDiaSearchParameters
    {
        // ── Fragment Extraction ──────────────────────────────────────────────

        /// <summary>Fragment m/z tolerance in ppm. Default: 20.</summary>
        public float PpmTolerance { get; set; } = 20f;

        /// <summary>
        /// Initial RT window half-width in minutes used for the calibration bootstrap pass.
        /// After calibration converges, adaptive windows (k×σ) are used automatically.
        /// Default: 5.0.
        /// </summary>
        public float RtToleranceMinutes { get; set; } = 5.0f;

        /// <summary>Minimum detected fragments for a result to be reported. Default: 3.</summary>
        public int MinFragmentsRequired { get; set; } = 3;

        /// <summary>Minimum score threshold (pre-FDR filter). Default: 0.0 (no pre-filtering).</summary>
        public float MinScoreThreshold { get; set; } = 0.0f;

        // ── Parallelism ──────────────────────────────────────────────────────

        /// <summary>Max extraction threads. -1 = all cores. Default: -1.</summary>
        public int MaxThreads { get; set; } = -1;

        /// <summary>Whether to attempt GPU acceleration. Default: false.</summary>
        public bool PreferGpu { get; set; } = false;

        // ── Calibration ──────────────────────────────────────────────────────

        /// <summary>
        /// Enable iterative RT calibration (DiaCalibrationPipeline).
        /// When true, produces substantially more IDs (benchmark: 29,157 vs ~16,000 without).
        /// Disable only for debugging or when the library has no meaningful RT values.
        /// Default: true.
        /// </summary>
        public bool UseIrtCalibration { get; set; } = true;

        // ── Output Control ───────────────────────────────────────────────────
        // These are MetaMorpheus-specific and are NOT passed to mzLib.

        /// <summary>Write per-file and summary diagnostic files. Default: true.</summary>
        public bool WriteDiagnostics { get; set; } = true;

        /// <summary>Include decoy results in output TSV files. Default: false.</summary>
        public bool WriteDecoyResults { get; set; } = false;

        /// <summary>Include contaminant results in output TSV files. Default: true.</summary>
        public bool WriteContaminantResults { get; set; } = true;

        /// <summary>Write per-file result files in addition to the combined results. Default: true.</summary>
        public bool WriteIndividualFiles { get; set; } = true;

        /// <summary>Include results above the q-value threshold in output. Default: true.</summary>
        public bool WriteHighQValueResults { get; set; } = true;
        public bool WritePsmTsv { get; set; } = true;
        /// <summary>
        /// Path to a spectral library file (.msp) to use for DIA searching.
        /// When empty, library spectra are expected to be generated from Koina/Prosit predictions.
        /// </summary>
        public string SpectralLibraryPath { get; set; } = "";

        // ── Factory ──────────────────────────────────────────────────────────

        /// <summary>
        /// Creates the mzLib DiaSearchParameters from this configuration.
        /// Only passes parameters that DiaSearchParameters actually exposes.
        /// Calibration settings (UseIrtCalibration) are handled at the DiaEngine level,
        /// not inside DiaSearchParameters.
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
                PreferGpu = PreferGpu
            };
        }
    }
}