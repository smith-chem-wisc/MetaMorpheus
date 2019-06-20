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

        /// <summary>
        /// At the time when the ~10% of the data gets chosen for training, another 10% gets chosen for evaluation. Then after training, the effectiveness of the 
        /// model gets evaluated on the test set. The results of that evaluation are converted to text values called BinarySearchTreeMetrics and this gets written to the results.tsv
        /// </summary>
        public string BinarySearchTreeMetrics { get; set; }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine($"{AnalysisType}s within 1% FDR: {PsmsWithin1PercentFdr.ToString()}");
            sb.AppendLine($"Delta Score Used for FDR Analysis: {DeltaScoreImprovement.ToString()}");
            sb.AppendLine(BinarySearchTreeMetrics);
            return sb.ToString();
        }
    }
}