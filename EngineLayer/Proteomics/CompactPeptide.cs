using Proteomics;
using System;
using System.Linq;
using System.Collections.Generic;

namespace EngineLayer
{
    [Serializable]
    public class CompactPeptide : CompactPeptideBase
    {

        #region Public Constructors

        public CompactPeptide(PeptideWithSetModifications peptideWithSetModifications)
        {
            double theMass = 0;
            ModificationWithMass pep_n_term_variable_mod;
            if (peptideWithSetModifications.allModsOneIsNterminus.TryGetValue(1, out pep_n_term_variable_mod))
                foreach (double nl in pep_n_term_variable_mod.neutralLosses)
                    theMass = pep_n_term_variable_mod.monoisotopicMass - nl;
            else
                theMass = 0;
            NTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, theMass, 1, 1).ToArray();

            ModificationWithMass pep_c_term_variable_mod;
            theMass = 0;
            if (peptideWithSetModifications.allModsOneIsNterminus.TryGetValue(peptideWithSetModifications.Length + 2, out pep_c_term_variable_mod))
                foreach (double nl in pep_c_term_variable_mod.neutralLosses)
                    theMass = pep_c_term_variable_mod.monoisotopicMass - nl;
            else
                theMass = 0;
            CTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, theMass, peptideWithSetModifications.Length, -1).ToArray();

            MonoisotopicMassIncludingFixedMods = peptideWithSetModifications.MonoisotopicMass;
        }

        #endregion Public Constructors

    }
}