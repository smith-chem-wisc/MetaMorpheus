using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.NonSpecificEnzymeSearch
{
    public class NonSpecificEnzymeSequencesToActualPeptidesResults : MetaMorpheusEngineResults
    {
        #region Public Constructors

        public NonSpecificEnzymeSequencesToActualPeptidesResults(MetaMorpheusEngine s, Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching) : base(s)
        {
            this.compactPeptideToProteinPeptideMatching = compactPeptideToProteinPeptideMatching;
            
        }

        #endregion Public Constructors

        #region Public Properties
        
        public Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching { get; }
        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            //i dont know if this is exactly what I want to do here
            sb.Append("CompactPeptideToProteinPeptideMatching.Count: " + compactPeptideToProteinPeptideMatching.Count);
            return sb.ToString();
        }

        #endregion Public Methods
    }
}

