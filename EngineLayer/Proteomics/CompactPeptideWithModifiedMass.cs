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

        #region Public Properties

        public override double[] CTerminalMasses { get; }
        public override double[] NTerminalMasses { get; }
        public override double MonoisotopicMassIncludingFixedMods { get; }

        #endregion Public Properties
    }
}
