using System;
using System.Linq;

namespace EngineLayer
{
    [Serializable]
    public class CompactPeptide : CompactPeptideBase
    {
        public CompactPeptide(PeptideWithSetModifications peptideWithSetModifications, TerminusType terminusType)
        {
            NTerminalMasses = null;
            CTerminalMasses = null;
            if (terminusType == TerminusType.None || terminusType == TerminusType.N)
            {
                NTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, 0, 0, 1).ToArray();
            }
            if (terminusType == TerminusType.None || terminusType == TerminusType.C)
            {
                CTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, 0, peptideWithSetModifications.Length + 1, -1).ToArray();
            }
            MonoisotopicMassIncludingFixedMods = peptideWithSetModifications.MonoisotopicMass;
        }
    }
}