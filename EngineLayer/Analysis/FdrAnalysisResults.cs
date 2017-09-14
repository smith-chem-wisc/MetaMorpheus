using System.Text;

namespace EngineLayer.Analysis
{
    public class FdrAnalysisResults : MetaMorpheusEngineResults
    {
        #region Public Constructors

        public FdrAnalysisResults(FdrAnalysisEngine s) : base(s)
        {
        }

        public int PsmsWithin1percentFdr { get; set; }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.Append("PSMs within 1% fdr: " + PsmsWithin1percentFdr);
            return sb.ToString();
        }

        #endregion Public Methods
    }
}