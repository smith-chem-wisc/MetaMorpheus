using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace IndexSearchAndAnalyze
{
    public class AnalysisResults : MyResults
    {
        public List<NewPsmWithFDR>[] allResultingIdentifications { get; private set; }

        public AnalysisResults(AnalysisParams s) : base(s)
        {
        }

        public AnalysisResults(AnalysisParams s, List<NewPsmWithFDR>[] yeah) : this(s)
        {
            this.allResultingIdentifications = yeah;
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("AnalysisResults: ");
            sb.AppendLine();
            sb.Append(base.ToString());
            sb.AppendLine();
            sb.Append("All PSMS within 1% FDR: " + string.Join(", ", allResultingIdentifications.Select(b => b.Where(c => c.QValue <= 0.01).Count())));
            sb.AppendLine();

            return sb.ToString();
        }
    }
}