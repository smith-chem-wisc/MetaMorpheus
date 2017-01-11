using InternalLogicEngineLayer;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InternalLogicEngineLayer
{
    public class AnalysisResults : MyResults
    {
        public List<NewPsmWithFDR>[] allResultingIdentifications { get; private set; }

        public AnalysisResults(AnalysisEngine s) : base(s)
        {
        }

        public AnalysisResults(AnalysisEngine s, List<NewPsmWithFDR>[] yeah) : this(s)
        {
            this.allResultingIdentifications = yeah;
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append(base.ToString());
            sb.AppendLine();
            sb.Append("All PSMS within 1% FDR: " + string.Join(", ", allResultingIdentifications.Select(b => b.Where(c => c.QValue <= 0.01).Count())));

            return sb.ToString();
        }
    }
}