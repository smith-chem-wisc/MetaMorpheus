using Proteomics;
using System;
using System.Linq;

namespace EngineLayer
{
    [Serializable]
    public class CompactPeptide : CompactPeptideBase
    {
        #region Public Constructors

        public CompactPeptide(PeptideWithSetModifications peptideWithSetModifications, TerminusType terminusType)
        {
            NTerminalMasses = null;
            CTerminalMasses = null;
            if (terminusType == TerminusType.None || terminusType == TerminusType.N)
            {
                double theMass = 0;
                if (peptideWithSetModifications.allModsOneIsNterminus.TryGetValue(1, out ModificationWithMass pep_n_term_variable_mod))
                    foreach (double nl in pep_n_term_variable_mod.neutralLosses)
                        theMass = pep_n_term_variable_mod.monoisotopicMass - nl;
                else
                    theMass = 0;
                NTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, theMass, 1, 1).ToArray();
            }
            if (terminusType == TerminusType.None || terminusType == TerminusType.C)
            {
                double theMass = 0;
                if (peptideWithSetModifications.allModsOneIsNterminus.TryGetValue(peptideWithSetModifications.Length + 2, out ModificationWithMass pep_c_term_variable_mod))
                    foreach (double nl in pep_c_term_variable_mod.neutralLosses)
                        theMass = pep_c_term_variable_mod.monoisotopicMass - nl;
                else
                    theMass = 0;
                CTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, theMass, peptideWithSetModifications.Length, -1).ToArray();
            }
            MonoisotopicMassIncludingFixedMods = peptideWithSetModifications.MonoisotopicMass;
        }

        #endregion Public Constructors
    }
}