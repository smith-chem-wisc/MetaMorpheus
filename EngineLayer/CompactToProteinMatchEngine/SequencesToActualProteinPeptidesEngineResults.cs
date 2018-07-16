using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer
{
    public class SequencesToActualProteinPeptidesEngineResults : MetaMorpheusEngineResults
    {
        public SequencesToActualProteinPeptidesEngineResults(MetaMorpheusEngine s, Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching) : base(s)
        {
            CompactPeptideToProteinPeptideMatching = compactPeptideToProteinPeptideMatching;
        }

        public Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> CompactPeptideToProteinPeptideMatching { get; }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.Append("CompactPeptideToProteinPeptideMatching.Count: " + CompactPeptideToProteinPeptideMatching.Count);
            return sb.ToString();
        }
    }
}