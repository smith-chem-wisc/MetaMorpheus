// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using System.Text;

namespace EngineLayer.FdrAnalysisDia
{
    /// <summary>
    /// Results summary from DIA FDR analysis.
    /// Analogous to FdrAnalysisResults for DDA, but tailored to DIA precursor matches
    /// which have no notch-level ambiguity or protease grouping.
    /// </summary>
    public class FdrAnalysisResultsDia : MetaMorpheusEngineResults
    {
        public FdrAnalysisResultsDia(FdrAnalysisEngineDia engine, string analysisType) : base(engine)
        {
            AnalysisType = analysisType;
        }

        /// <summary>Number of target results passing 1% global q-value threshold</summary>
        public int ResultsWithin1PercentFdr { get; set; }

        /// <summary>Number of unique peptide sequences passing 1% peptide-level q-value threshold</summary>
        public int PeptidesWithin1PercentFdr { get; set; }

        /// <summary>Total target results (before FDR filtering)</summary>
        public int TotalTargets { get; set; }

        /// <summary>Total decoy results (before FDR filtering)</summary>
        public int TotalDecoys { get; set; }

        /// <summary>Description of the analysis level</summary>
        private string AnalysisType { get; set; }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine($"DIA {AnalysisType} FDR Analysis:");
            sb.AppendLine($"  Total targets: {TotalTargets}");
            sb.AppendLine($"  Total decoys: {TotalDecoys}");
            sb.AppendLine($"  {AnalysisType}s within 1% FDR: {ResultsWithin1PercentFdr}");
            sb.AppendLine($"  Peptides within 1% FDR: {PeptidesWithin1PercentFdr}");
            return sb.ToString();
        }
    }
}
