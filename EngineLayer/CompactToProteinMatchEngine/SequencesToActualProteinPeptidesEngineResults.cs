using System.Collections.Generic;
using System.Text;

namespace EngineLayer
{
    public class SequencesToActualProteinPeptidesEngineResults : MetaMorpheusEngineResults
    {
        
        #region Public Constructors

        public SequencesToActualProteinPeptidesEngineResults(MetaMorpheusEngine s, Dictionary<Protease, Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>> proteaseSortedCompactPeptideToProteinPeptideMatching) : base(s)
        {
            this.proteaseSortedCompactPeptideToProteinPeptideMatching = proteaseSortedCompactPeptideToProteinPeptideMatching;
            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            foreach (var proteaseSet in proteaseSortedCompactPeptideToProteinPeptideMatching)
            {
                var CPWM = proteaseSet.Value;
                foreach (var CPWMkvp in CPWM)
                {
                    compactPeptideToProteinPeptideMatching.Add(CPWMkvp.Key, CPWMkvp.Value);
                }
            }
            this.compactPeptideToProteinPeptideMatching = compactPeptideToProteinPeptideMatching;
        }

        #endregion Public Constructors

        #region Public Properties
        public Dictionary<Protease, Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>> proteaseSortedCompactPeptideToProteinPeptideMatching { get; }
        public Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching { get; }
        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            //i dont know if this is exactly what I want to do here
            sb.Append("CompactPeptideToProteinPeptideMatching.Count: " + proteaseSortedCompactPeptideToProteinPeptideMatching.Count);
            return sb.ToString();
        }

        #endregion Public Methods
    }
}