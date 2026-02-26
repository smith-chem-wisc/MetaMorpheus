// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using EngineLayer;
using EngineLayer.FdrAnalysisDia;
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
        /// Which score type to use for FDR ranking.
        /// Default: SpectralAngle (generally better discriminative power than dot product).
        /// </summary>
        public DiaFdrScoreType FdrScoreType { get; set; } = DiaFdrScoreType.SpectralAngle;

        /// <summary>
        /// Q-value threshold for reporting "passing" results.
        /// Default: 0.01 (1% FDR).
        /// </summary>
        public double QValueThreshold { get; set; } = 0.01;

        /// <summary>Whether to write individual per-file result files</summary>
        public bool WriteIndividualFiles { get; set; } = true;

        /// <summary>Whether to include decoy results in output files</summary>
        public bool WriteDecoys { get; set; } = false;

        /// <summary>Whether to include results above the q-value threshold in output</summary>
        public bool WriteHighQValueResults { get; set; } = false;
    }
}