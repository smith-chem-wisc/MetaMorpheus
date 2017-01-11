using InternalLogicEngineLayer;
using OldInternalLogic;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InternalLogicEngineLayer
{
    public class AnalysisResults : MyResults
    {
        public Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> dict { get; private set; }

        public List<NewPsmWithFDR>[] allResultingIdentifications { get; private set; }

        public AnalysisResults(AnalysisEngine s, List<NewPsmWithFDR>[] yeah, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> dict) : base(s)
        {
            this.allResultingIdentifications = yeah;
            this.dict = dict;
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