using System.Collections.Generic;
using EngineLayer;

namespace EngineLayer
{
    public class SequencesToActualProteinPeptidesEngineResults : MetaMorpheusEngineResults
    {
        public Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> CompactPeptideToProteinPeptideMatching { get; }
        
        public SequencesToActualProteinPeptidesEngineResults(MetaMorpheusEngine s, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching) : base(s)
        {
            this.CompactPeptideToProteinPeptideMatching = compactPeptideToProteinPeptideMatching;
        }
    }
}