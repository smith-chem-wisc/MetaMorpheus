using System.Text;

namespace IndexSearchAndAnalyze
{
    public class AnalysisResults : MyResults
    {
        public AnalysisResults(AnalysisParams s) : base(s)
        {
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("AnalysisResults: ");
            sb.Append(base.ToString());
            return sb.ToString();
        }
    }
}