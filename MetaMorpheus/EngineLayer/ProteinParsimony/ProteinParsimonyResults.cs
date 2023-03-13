using System.Collections.Generic;
using System.Text;

namespace EngineLayer
{
    public class ProteinParsimonyResults : MetaMorpheusEngineResults
    {
        public ProteinParsimonyResults(ProteinParsimonyEngine proteinAnalysisEngine) : base(proteinAnalysisEngine)
        {
        }

        public List<ProteinGroup> ProteinGroups { get; set; }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            return sb.ToString();
        }
    }
}