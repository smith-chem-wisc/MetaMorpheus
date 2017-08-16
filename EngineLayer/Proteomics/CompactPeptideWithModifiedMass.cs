using System;
using System.Collections.Generic;

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
            this.ModifiedMass = MonoisotopicMassIncludingFixedMods;
        }

        #endregion Public Constructors

        #region Public Properties

        public double ModifiedMass { get; set; }

        #endregion Public Properties

        #region Public Methods

        public void SwapMonoisotopicMassWithModifiedMass()
        {
            double tempDouble = this.MonoisotopicMassIncludingFixedMods;
            this.MonoisotopicMassIncludingFixedMods = this.ModifiedMass;
            this.ModifiedMass = tempDouble;
        }

        public void CropTerminalMasses(TerminusType terminusType)
        {
            List<double> tempList = new List<double>();
            if (terminusType == TerminusType.N)
            {
                for (int i = 0; i < NTerminalMasses.Length; i++)
                {
                    if (NTerminalMasses[i] < MonoisotopicMassIncludingFixedMods)
                    {
                        tempList.Add(NTerminalMasses[i]);
                    }
                    else
                    {
                        NTerminalMasses = tempList.ToArray();
                        break;
                    }
                }
            }
            else
            {
                for (int i = 0; i < CTerminalMasses.Length; i++)
                {
                    if (CTerminalMasses[i] < MonoisotopicMassIncludingFixedMods)
                    {
                        tempList.Add(CTerminalMasses[i]);
                    }
                    else
                    {
                        CTerminalMasses = tempList.ToArray();
                        break;
                    }
                }
            }
        }

        #endregion Public Methods
    }
}