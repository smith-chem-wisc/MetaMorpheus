using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public class ProteinScoringAndFdrResults : MetaMorpheusEngineResults
    {
        #region Public Fields

        public List<ProteinGroup> sortedAndScoredProteinGroups; //used outside the function, better keep public

        #endregion Public Fields

        #region Public Constructors

        public ProteinScoringAndFdrResults(ProteinScoringAndFdrEngine proteinAnalysisEngine) : base(proteinAnalysisEngine)
        {
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.Append("Number of proteins within 1% FDR: " + sortedAndScoredProteinGroups.Count(b => b.QValue < 0.01));
            return sb.ToString();
        }

        #endregion Public Methods
    }
}