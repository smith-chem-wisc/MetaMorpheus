using System.Collections.Generic;
using System.Text;

namespace EngineLayer
{
    public class ProteinParsimonyResults : MetaMorpheusEngineResults
    {
        public int HypothesesAdded { get; set; }
        public int HypothesesRemoved { get; set; }
        public ProteinParsimonyResults(ProteinParsimonyEngine proteinAnalysisEngine) : base(proteinAnalysisEngine)
        {
        }

        public List<ProteinGroup> ProteinGroups { get; set; }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine($"Hypotheses Added to Spectral Matches:       {HypothesesAdded}");
            sb.AppendLine($"Hypotheses Removed from Spectral Matches:   {HypothesesRemoved}");
            return sb.ToString();
        }
    }
}
