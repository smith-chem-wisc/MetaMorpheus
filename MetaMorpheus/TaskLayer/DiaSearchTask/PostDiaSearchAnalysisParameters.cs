// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using EngineLayer;
using MassSpectrometry.Dia;
using System.Collections.Generic;

namespace TaskLayer
{
    /// <summary>
    /// Parameters container for DIA post-search analysis.
    /// 
    /// Analogous to PostSearchAnalysisParameters for DDA, but carries DIA-specific data:
    /// - DiaSearchResult list instead of SpectralMatch list
    /// - DiaSearchParameters instead of SearchParameters
    /// - No protease/modification context (DIA uses library-based identification)
    /// - No FlashLFQ results (DIA quantification is handled differently)
    /// </summary>
    public class PostDiaSearchAnalysisParameters
    {
        /// <summary>The task results object to populate with summary information</summary>
        public MyTaskResults SearchTaskResults { get; set; }

        /// <summary>Task ID for status reporting</summary>
        public string SearchTaskId { get; set; }

        /// <summary>DIA-specific search parameters used during extraction</summary>
        public DiaSearchParameters DiaSearchParameters { get; set; }

        /// <summary>All DIA precursor matches from the search (targets + decoys, pre-FDR)</summary>
        public List<DiaSearchResult> AllDiaResults { get; set; }

        /// <summary>Output folder for results files</summary>
        public string OutputFolder { get; set; }

        /// <summary>Output subfolder for per-file results</summary>
        public string IndividualResultsOutputFolder { get; set; }

        /// <summary>List of raw data files that were searched</summary>
        public List<string> CurrentRawFileList { get; set; }

        /// <summary>
        /// Which classifier to use for iterative FDR.
        /// Default: NeuralNetwork (production default, 29,157 IDs at 1% FDR on benchmark).
        /// LinearDiscriminant is the fast fallback (slightly fewer IDs, ~4x faster).
        /// </summary>
        public DiaClassifierType ClassifierType { get; set; } = DiaClassifierType.NeuralNetwork;

        /// <summary>
        /// Q-value threshold for reporting "passing" results.
        /// Default: 0.01 (1% FDR).
        /// </summary>
        public double QValueThreshold { get; set; } = 0.01;

        /// <summary>Whether to write individual per-file result files</summary>
        public bool WriteIndividualFiles { get; set; } = true;

        /// <summary>Whether to include decoy results in output files</summary>
        public bool WriteDecoys { get; set; } = false;
        /// <summary>
        /// Whether to write .psmtsv output files compatible with MetaDraw and
        /// downstream MetaMorpheus tools. Default: true.
        /// </summary>
        public bool WritePsmTsv { get; set; } = true;
        /// <summary>Whether to include results above the q-value threshold in output</summary>
        public bool WriteHighQValueResults { get; set; } = false;

        /// <summary>
        /// The scan index built during DiaEngine execution.
        /// Required by DiaFeatureExtractor for MS1 feature computation (isotope scores,
        /// MS1-MS2 correlation, precursor XIC intensity, precursor elution score).
        /// Set by DiaSearchTask after DiaEngine.Run() returns — the index must remain
        /// alive until PostDiaSearchAnalysisTask.Run() completes.
        /// Null if MS1 features are not available (feature computation will skip MS1 features).
        /// </summary>
        public DiaScanIndex DiaScanIndex { get; set; } = null;
    }
}