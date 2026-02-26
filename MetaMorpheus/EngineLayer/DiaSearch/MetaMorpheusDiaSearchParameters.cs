// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using MassSpectrometry.Dia;

namespace EngineLayer.DiaSearch
{
    /// <summary>
    /// TOML-serializable parameter class for the MetaMorpheus DIA search task layer.
    /// Wraps mzLib's <see cref="DiaSearchParameters"/> with additional task-level settings
    /// (spectral library path, output options, etc.).
    /// 
    /// Nett serializes/deserializes this class automatically for .toml task files.
    /// Only TOML-friendly types are used (primitives, strings).
    /// </summary>
    public class MetaMorpheusDiaSearchParameters
    {
        #region Fragment Extraction (maps to mzLib DiaSearchParameters)

        /// <summary>
        /// Fragment m/z tolerance in parts-per-million for XIC extraction.
        /// Default: 20 ppm (typical for Orbitrap MS2).
        /// </summary>
        public float PpmTolerance { get; set; } = 20f;

        /// <summary>
        /// Half-width of the RT window around the predicted/library retention time (in minutes).
        /// Default: 5.0 minutes.
        /// </summary>
        public float RtToleranceMinutes { get; set; } = 5.0f;

        /// <summary>
        /// Minimum number of fragment ions that must yield XIC data points
        /// for a precursor to be reported as a match.
        /// Default: 3.
        /// </summary>
        public int MinFragmentsRequired { get; set; } = 3;

        /// <summary>
        /// Minimum dot product score threshold. Precursors below this are excluded.
        /// Default: 0.0 (no filtering; leave to downstream FDR).
        /// </summary>
        public float MinScoreThreshold { get; set; } = 0.0f;

        /// <summary>
        /// Whether to prefer GPU extraction when a compatible device is detected.
        /// Falls back to CPU if no GPU is available.
        /// Default: false.
        /// </summary>
        public bool PreferGpu { get; set; } = false;

        /// <summary>
        /// Maximum threads for window-level parallelism.
        /// -1 = use Environment.ProcessorCount.
        /// Default: -1.
        /// </summary>
        public int MaxThreads { get; set; } = -1;

        #endregion

        #region Task-Level Settings

        /// <summary>
        /// Path to the spectral library file (.msp format).
        /// If empty or null, the task will look for a library in the database file list.
        /// </summary>
        public string SpectralLibraryPath { get; set; } = "";

        /// <summary>
        /// Whether to include decoy identifications in the output TSV.
        /// Default: false (targets only).
        /// </summary>
        public bool WriteDecoyResults { get; set; } = false;

        /// <summary>
        /// Whether to write a diagnostics file with timing and throughput information.
        /// Default: true.
        /// </summary>
        public bool WriteDiagnostics { get; set; } = true;

        #endregion

        /// <summary>
        /// Converts task-level parameters to the mzLib engine-level parameters.
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
