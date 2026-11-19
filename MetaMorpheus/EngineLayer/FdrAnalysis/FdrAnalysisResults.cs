using System.Text;

namespace EngineLayer.FdrAnalysis
{
    public class FdrAnalysisResults : MetaMorpheusEngineResults
    {
        public FdrAnalysisResults(FdrAnalysisEngine s, string analysisType) : base(s)
        {
            DeltaScoreImprovement = false;
            AnalysisType = analysisType;
        }

        /// <summary>
        /// The target matches that pass the search's identification filter at its threshold (<see cref="IdentificationFilter"/>),
        /// 1% by default.
        /// </summary>
        public int PsmsWithin1PercentFdr { get; set; }

        /// <summary>The filter and threshold <see cref="PsmsWithin1PercentFdr"/> counts by: "q-value notch &lt;= 0.01".</summary>
        public string IdentificationFilter { get; set; }
        public bool DeltaScoreImprovement { get; set; }
        private string AnalysisType { get; set; }
        public string BinarySearchTreeMetrics { get; set; } //See PEPValueAnalysisGeneric public static string PrintBinaryClassificationMetrics method for explanation

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine($"{AnalysisType}s with {IdentificationFilter}: {PsmsWithin1PercentFdr}");
            sb.AppendLine($"Delta Score Used for FDR Analysis: {DeltaScoreImprovement.ToString()}");
            sb.AppendLine(BinarySearchTreeMetrics);
            return sb.ToString();
        }
    }
}