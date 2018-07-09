using System.Text;

namespace EngineLayer.FdrAnalysis
{
    public class FdrAnalysisResults : MetaMorpheusEngineResults
    {
        public FdrAnalysisResults(FdrAnalysisEngine s) : base(s)
        {
            DeltaScoreImprovement = false;
        }

        public int PsmsWithin1PercentFdr { get; set; }
        public bool DeltaScoreImprovement { get; set; }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("PSMs within 1% fdr: " + PsmsWithin1PercentFdr);
            sb.AppendLine("Delta Score Used for FDR Analysis: " + DeltaScoreImprovement);
            return sb.ToString();
        }
    }
}