using MetaMorpheus;
using System.Collections.Generic;

namespace GoodGUI
{
    internal class SearchParamsObject
    {
        public List<double> accepted_precursor_mass_errors;
        public InitiatorMethionineBehavior initiatorMethionineBehavior;
        public MassTolerance precursorMassTolerance;
        public MassTolerance productMassTolerance;
        public Protease protease;
        public bool lookAtDecoys;
        public int missedCleavages;
        public int max_variable_mod_isoforms;
        public bool bions;
        public bool yions;
        public int max_mods_for_peptide;

        public SearchParamsObject(bool lookAtDecoys, Protease protease, int missedCleavages, InitiatorMethionineBehavior initiatorMethionineBehavior, MassTolerance precursorMassTolerance, List<double> accepted_precursor_mass_errors, MassTolerance productMassTolerance, bool bions, bool yions, int max_variable_mod_isoforms, int max_mods_for_peptide)
        {
            this.lookAtDecoys = lookAtDecoys;
            this.protease = protease;
            this.missedCleavages = missedCleavages;
            this.initiatorMethionineBehavior = initiatorMethionineBehavior;
            this.precursorMassTolerance = precursorMassTolerance;
            this.accepted_precursor_mass_errors = accepted_precursor_mass_errors;
            this.productMassTolerance = productMassTolerance;
            this.bions = bions;
            this.yions = yions;
            this.max_variable_mod_isoforms = max_variable_mod_isoforms;
            this.max_mods_for_peptide = max_mods_for_peptide;
        }
    }
}