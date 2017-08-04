using System;

namespace EngineLayer
{
    [Serializable]
    internal class CompactPeptideWithModifiedMass : CompactPeptideBase
    {

        #region Public Constructors

        public CompactPeptideWithModifiedMass(CompactPeptideBase cp, double MonoisotopicMassIncludingFixedMods)
        {
            this.CTerminalMasses = cp.CTerminalMasses;
            this.NTerminalMasses = cp.NTerminalMasses;
            this.MonoisotopicMassIncludingFixedMods = cp.MonoisotopicMassIncludingFixedMods;
            this.modifiedMass = MonoisotopicMassIncludingFixedMods;
        }

        #endregion Public Constructors

        #region Public Properties

        public double modifiedMass { get; set; }

        #endregion Public Properties

        #region Public Methods

        public void AssignCorrectMass()
        {
            this.MonoisotopicMassIncludingFixedMods = this.modifiedMass;
        }

        #endregion Public Methods

    }
}