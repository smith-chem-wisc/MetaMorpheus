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

        public int PsmsWithin1PercentFdr { get; set; }
        public bool DeltaScoreImprovement { get; set; }
        private string AnalysisType { get; set; }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine($"{AnalysisType}s within 1% FDR: {PsmsWithin1PercentFdr.ToString()}");
            sb.AppendLine($"Delta Score Used for FDR Analysis: {DeltaScoreImprovement.ToString()}");
            return sb.ToString();
        }
    }
}