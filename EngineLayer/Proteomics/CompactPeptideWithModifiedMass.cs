using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer
{
    [Serializable]
    class CompactPeptideWithModifiedMass : CompactPeptideBase
    {
        #region Public Constructors
        public double modifiedMass { get; set; }

        public CompactPeptideWithModifiedMass(CompactPeptideBase cp, double MonoisotopicMassIncludingFixedMods)
        {
            this.CTerminalMasses = cp.CTerminalMasses;
            this.NTerminalMasses = cp.NTerminalMasses;
            this.MonoisotopicMassIncludingFixedMods = cp.MonoisotopicMassIncludingFixedMods;
            this.modifiedMass = MonoisotopicMassIncludingFixedMods;
        }

        #endregion Public Constructors

        public void AssignCorrectMass()
        {
            this.MonoisotopicMassIncludingFixedMods = this.modifiedMass;
        }
    }
}
