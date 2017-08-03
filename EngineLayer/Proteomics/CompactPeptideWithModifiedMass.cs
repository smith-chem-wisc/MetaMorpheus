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

        public CompactPeptideWithModifiedMass(CompactPeptide cp, double MonoisotopicMassIncludingFixedMods)
        {
            this.CTerminalMasses = cp.CTerminalMasses;
            this.NTerminalMasses = cp.NTerminalMasses;
            this.MonoisotopicMassIncludingFixedMods = MonoisotopicMassIncludingFixedMods;
        }

        #endregion Public Constructors
    }
}
