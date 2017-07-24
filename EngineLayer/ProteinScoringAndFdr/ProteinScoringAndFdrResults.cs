using System.Collections.Generic;

namespace EngineLayer
{
    public class ProteinScoringAndFdrResults : MetaMorpheusEngineResults
    {

        #region Public Fields

        public List<ProteinGroup> sortedAndScoredProteinGroups;

        #endregion Public Fields

        #region Public Constructors

        public ProteinScoringAndFdrResults(ProteinScoringAndFdrEngine proteinAnalysisEngine) : base(proteinAnalysisEngine)
        {
        }

        #endregion Public Constructors

    }
}