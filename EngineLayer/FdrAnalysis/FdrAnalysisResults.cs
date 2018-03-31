using System.Text;

namespace EngineLayer.FdrAnalysis
{
    public class FdrAnalysisResults : MetaMorpheusEngineResults
    {
        #region Public Constructors

        public FdrAnalysisResults(FdrAnalysisEngine s) : base(s)
        {
            DeltaScoreImprovement = false;
        }

        #endregion Public Constructors

        #region Public Properties

        public int PsmsWithin1PercentFdr { get; set; }
        public bool DeltaScoreImprovement { get; set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.Append("PSMs within 1% fdr: " + PsmsWithin1PercentFdr);
            sb.Append("Delta Score Improvement: " + DeltaScoreImprovement);
            return sb.ToString();
        }

        #endregion Public Methods
    }
}