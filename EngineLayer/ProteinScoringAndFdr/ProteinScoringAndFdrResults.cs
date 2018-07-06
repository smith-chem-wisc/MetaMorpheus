using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public class ProteinScoringAndFdrResults : MetaMorpheusEngineResults
    {
        public List<ProteinGroup> SortedAndScoredProteinGroups;

        public ProteinScoringAndFdrResults(ProteinScoringAndFdrEngine proteinAnalysisEngine) : base(proteinAnalysisEngine)
        {
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.Append("Number of proteins within 1% FDR: " + SortedAndScoredProteinGroups.Count(b => b.QValue < 0.01));
            return sb.ToString();
        }
    }
}