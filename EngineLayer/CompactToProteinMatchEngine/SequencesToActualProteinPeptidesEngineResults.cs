using System.Collections.Generic;
using System.Text;

namespace EngineLayer
{
    public class SequencesToActualProteinPeptidesEngineResults : MetaMorpheusEngineResults
    {

        #region Public Constructors

        public SequencesToActualProteinPeptidesEngineResults(MetaMorpheusEngine s, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching) : base(s)
        {
            this.CompactPeptideToProteinPeptideMatching = compactPeptideToProteinPeptideMatching;
        }

        #endregion Public Constructors

        #region Public Properties

        public Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> CompactPeptideToProteinPeptideMatching { get; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.Append("CompactPeptideToProteinPeptideMatching.Count: " + CompactPeptideToProteinPeptideMatching.Count);
            return sb.ToString();
        }

        #endregion Public Methods

    }
}