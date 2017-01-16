using System;
using System.Collections.Generic;
using System.Linq;

namespace OldInternalLogic
{
    public class PeptideWithPossibleModifications : Peptide
    {
        #region Public Fields

        public Dictionary<int, List<MorpheusModification>> twoBasedFixedModificationss = new Dictionary<int, List<MorpheusModification>>();

        #endregion Public Fields

        #region Public Constructors

        public PeptideWithPossibleModifications(int OneBasedStartResidueNumberInProtein, int OneBasedEndResidueNumberInProtein, Protein ParentProtein, int missedCleavages, string PeptideDescription)
            : base(ParentProtein, OneBasedStartResidueNumberInProtein, OneBasedEndResidueNumberInProtein)
        {
            OneBasedPossibleLocalizedModifications = protein.OneBasedPossibleLocalizedModifications.Where(ok => ok.Key >= OneBasedStartResidueInProtein && ok.Key <= OneBasedEndResidueInProtein).ToDictionary(ok => ok.Key - OneBasedStartResidueInProtein + 1, ok => ok.Value);
            this.missedCleavages = missedCleavages;
            this.PeptideDescription = PeptideDescription;
        }

        #endregion Public Constructors

        #region Public Properties

        public Dictionary<int, List<MorpheusModification>> OneBasedPossibleLocalizedModifications { get; private set; }
        public int missedCleavages { get; private set; }

        public Dictionary<int, List<MorpheusModification>> FixedModifications { get; private set; }

        public Dictionary<int, MorpheusModification> VariableModifications { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public void SetFixedModifications(IEnumerable<MorpheusModification> allKnownFixedModifications)
        {
            twoBasedFixedModificationss = new Dictionary<int, List<MorpheusModification>>();
            for (int i = 0; i <= Length + 3; i++)
            {
                foreach (MorpheusModification mod in allKnownFixedModifications)
                {
                    if (i == 0 && (OneBasedStartResidueInProtein == 1 || OneBasedStartResidueInProtein == 2))
                    {
                        if (mod.Type == ModificationType.ProteinNTerminus && (mod.AminoAcid.Equals(this[0]) || mod.AminoAcid.Equals('\0')))
                        {
                            List<MorpheusModification> val;
                            if (twoBasedFixedModificationss.TryGetValue(i, out val))
                                val.Add(mod);
                            else
                                twoBasedFixedModificationss.Add(i, new List<MorpheusModification> { mod });
                        }
                    }
                    else if (i == 1)
                    {
                        if (mod.Type == ModificationType.PeptideNTerminus && (mod.AminoAcid.Equals(this[0]) || mod.AminoAcid.Equals('\0')))
                        {
                            List<MorpheusModification> val;
                            if (twoBasedFixedModificationss.TryGetValue(i, out val))
                                val.Add(mod);
                            else
                                twoBasedFixedModificationss.Add(i, new List<MorpheusModification> { mod });
                        }
                    }
                    else if (i >= 2 && i <= Length + 1)
                    {
                        if (mod.Type == ModificationType.AminoAcidResidue && (mod.AminoAcid.Equals(this[i - 2]) || mod.AminoAcid.Equals('\0')))
                        {
                            List<MorpheusModification> val;
                            if (twoBasedFixedModificationss.TryGetValue(i, out val))
                                val.Add(mod);
                            else
                                twoBasedFixedModificationss.Add(i, new List<MorpheusModification> { mod });
                        }
                    }
                    else if (i == Length + 2)
                    {
                        if (mod.Type == ModificationType.PeptideCTerminus && (mod.AminoAcid.Equals(this[Length]) || mod.AminoAcid.Equals('\0')))
                        {
                            List<MorpheusModification> val;
                            if (twoBasedFixedModificationss.TryGetValue(i, out val))
                                val.Add(mod);
                            else
                                twoBasedFixedModificationss.Add(i, new List<MorpheusModification> { mod });
                        }
                    }
                    else if (i == Length + 3 && OneBasedEndResidueInProtein == protein.Length)
                    {
                        if (mod.Type == ModificationType.ProteinCTerminus && (mod.AminoAcid.Equals(this[Length]) || mod.AminoAcid.Equals('\0')))
                        {
                            List<MorpheusModification> val;
                            if (twoBasedFixedModificationss.TryGetValue(i, out val))
                                val.Add(mod);
                            else
                                twoBasedFixedModificationss.Add(i, new List<MorpheusModification> { mod });
                        }
                    }
                }
            }
        }

        public IEnumerable<PeptideWithSetModifications> GetPeptideWithSetModifications(List<MorpheusModification> variableModifications, int maximumVariableModificationIsoforms, int max_mods_for_peptide)
        {
            var two_based_possible_variable_and_localizeable_modifications = new Dictionary<int, UniqueModificationsList>(Length + 4);

            var prot_n_term_variable_mods = new UniqueModificationsList();
            two_based_possible_variable_and_localizeable_modifications.Add(0, prot_n_term_variable_mods);

            var pep_n_term_variable_mods = new UniqueModificationsList();
            two_based_possible_variable_and_localizeable_modifications.Add(1, pep_n_term_variable_mods);

            var pep_c_term_variable_mods = new UniqueModificationsList();
            two_based_possible_variable_and_localizeable_modifications.Add(Length + 2, pep_c_term_variable_mods);

            var prot_c_term_variable_mods = new UniqueModificationsList();
            two_based_possible_variable_and_localizeable_modifications.Add(Length + 3, prot_c_term_variable_mods);

            foreach (MorpheusModification variable_modification in variableModifications)
            {
                if (variable_modification.Type == ModificationType.ProteinNTerminus && (OneBasedStartResidueInProtein == 1 || (OneBasedStartResidueInProtein == 2 && protein[0] == 'M'))
                    && (variable_modification.AminoAcid == char.MinValue || this[0] == variable_modification.AminoAcid))
                {
                    prot_n_term_variable_mods.Add(variable_modification);
                }

                if (variable_modification.Type == ModificationType.PeptideNTerminus && (variable_modification.AminoAcid == char.MinValue || this[0] == variable_modification.AminoAcid))
                {
                    pep_n_term_variable_mods.Add(variable_modification);
                }

                for (int r = 0; r < Length; r++)
                {
                    if (variable_modification.Type == ModificationType.AminoAcidResidue && this[r] == variable_modification.AminoAcid)
                    {
                        UniqueModificationsList residue_variable_mods;
                        if (!two_based_possible_variable_and_localizeable_modifications.TryGetValue(r + 2, out residue_variable_mods))
                        {
                            residue_variable_mods = new UniqueModificationsList();
                            residue_variable_mods.Add(variable_modification);
                            two_based_possible_variable_and_localizeable_modifications.Add(r + 2, residue_variable_mods);
                        }
                        else
                        {
                            residue_variable_mods.Add(variable_modification);
                        }
                    }
                }

                if (variable_modification.Type == ModificationType.PeptideCTerminus && (variable_modification.AminoAcid == char.MinValue || this[Length - 1] == variable_modification.AminoAcid))
                {
                    pep_c_term_variable_mods.Add(variable_modification);
                }

                if (variable_modification.Type == ModificationType.ProteinCTerminus && (OneBasedEndResidueInProtein == protein.Length - 1)
                    && (variable_modification.AminoAcid == char.MinValue || this[Length - 1] == variable_modification.AminoAcid))
                {
                    prot_c_term_variable_mods.Add(variable_modification);
                }
            }

            if (OneBasedPossibleLocalizedModifications != null)
            {
                foreach (KeyValuePair<int, List<MorpheusModification>> kvp in OneBasedPossibleLocalizedModifications)
                {
                    foreach (MorpheusModification variable_modification in kvp.Value)
                    {
                        if (((variable_modification.Type == ModificationType.ProteinNTerminus && !protein.isDecoy) ||
                            (variable_modification.Type == ModificationType.ProteinCTerminus && protein.isDecoy))
                            && (OneBasedStartResidueInProtein == 1 || (OneBasedStartResidueInProtein == 2 && protein[0] == 'M'))
                            && kvp.Key == 1)
                        {
                            prot_n_term_variable_mods.Add(variable_modification);
                        }

                        if (((variable_modification.Type == ModificationType.PeptideNTerminus && !protein.isDecoy) ||
                            (variable_modification.Type == ModificationType.PeptideCTerminus && protein.isDecoy))
                            && kvp.Key == 1)
                        {
                            pep_n_term_variable_mods.Add(variable_modification);
                        }

                        for (int r = 1; r <= Length; r++)
                        {
                            if (variable_modification.Type == ModificationType.AminoAcidResidue
                                && kvp.Key == r)
                            {
                                UniqueModificationsList residue_variable_mods;
                                if (!two_based_possible_variable_and_localizeable_modifications.TryGetValue(r + 1, out residue_variable_mods))
                                {
                                    residue_variable_mods = new UniqueModificationsList();
                                    residue_variable_mods.Add(variable_modification);
                                    two_based_possible_variable_and_localizeable_modifications.Add(r + 1, residue_variable_mods);
                                }
                                else
                                {
                                    residue_variable_mods.Add(variable_modification);
                                }
                            }
                        }

                        if (((variable_modification.Type == ModificationType.PeptideCTerminus && !protein.isDecoy) ||
                            (variable_modification.Type == ModificationType.PeptideNTerminus && protein.isDecoy))
                            && kvp.Key == Length)
                        {
                            pep_c_term_variable_mods.Add(variable_modification);
                        }

                        if (((variable_modification.Type == ModificationType.ProteinCTerminus && !protein.isDecoy) ||
                            (variable_modification.Type == ModificationType.ProteinNTerminus && protein.isDecoy))
                            && (OneBasedEndResidueInProtein == protein.Length)
                            && kvp.Key == Length)
                        {
                            prot_c_term_variable_mods.Add(variable_modification);
                        }
                    }
                }
            }

            int variable_modification_isoforms = 0;
            foreach (Dictionary<int, MorpheusModification> kvp in GetVariableModificationPatterns(two_based_possible_variable_and_localizeable_modifications, max_mods_for_peptide))
            {
                yield return new PeptideWithSetModifications(this, kvp);
                variable_modification_isoforms++;
                if (variable_modification_isoforms == maximumVariableModificationIsoforms)
                    yield break;
            }
        }

        #endregion Public Methods

        #region Protected Methods

        protected IEnumerable<Dictionary<int, MorpheusModification>> GetVariableModificationPatterns(Dictionary<int, UniqueModificationsList> possibleVariableModifications, int max_mods_for_peptide)
        {
            if (possibleVariableModifications.Count == 0)
            {
                yield return null;
            }
            else
            {
                var possible_variable_modifications = new Dictionary<int, UniqueModificationsList>(possibleVariableModifications);

                int[] base_variable_modification_pattern = new int[Length + 4];
                var totalAvailableMods = possible_variable_modifications.Select(b => b.Value == null ? 0 : b.Value.Count).Sum();
                for (int variable_modifications = 0; variable_modifications <= Math.Min(totalAvailableMods, max_mods_for_peptide); variable_modifications++)
                {
                    foreach (int[] variable_modification_pattern in GetVariableModificationPatterns(new List<KeyValuePair<int, UniqueModificationsList>>(possible_variable_modifications), possible_variable_modifications.Count - variable_modifications, base_variable_modification_pattern, 0))
                    {
                        yield return GetNewVariableModificationPattern(variable_modification_pattern, possible_variable_modifications);
                    }
                }
            }
        }

        #endregion Protected Methods

        #region Private Methods

        private static IEnumerable<int[]> GetVariableModificationPatterns(List<KeyValuePair<int, UniqueModificationsList>> possibleVariableModifications, int unmodifiedResiduesDesired, int[] variableModificationPattern, int index)
        {
            if (index < possibleVariableModifications.Count - 1)
            {
                if (unmodifiedResiduesDesired > 0)
                {
                    variableModificationPattern[possibleVariableModifications[index].Key] = 0;
                    foreach (int[] new_variable_modification_pattern in GetVariableModificationPatterns(possibleVariableModifications, unmodifiedResiduesDesired - 1, variableModificationPattern, index + 1))
                    {
                        yield return new_variable_modification_pattern;
                    }
                }
                if (unmodifiedResiduesDesired < possibleVariableModifications.Count - index)
                {
                    for (int i = 1; i <= possibleVariableModifications[index].Value.Count; i++)
                    {
                        variableModificationPattern[possibleVariableModifications[index].Key] = i;
                        foreach (int[] new_variable_modification_pattern in GetVariableModificationPatterns(possibleVariableModifications, unmodifiedResiduesDesired, variableModificationPattern, index + 1))
                        {
                            yield return new_variable_modification_pattern;
                        }
                    }
                }
            }
            else
            {
                if (unmodifiedResiduesDesired > 0)
                {
                    variableModificationPattern[possibleVariableModifications[index].Key] = 0;
                    yield return variableModificationPattern;
                }
                else
                {
                    for (int i = 1; i <= possibleVariableModifications[index].Value.Count; i++)
                    {
                        variableModificationPattern[possibleVariableModifications[index].Key] = i;
                        yield return variableModificationPattern;
                    }
                }
            }
        }

        private static Dictionary<int, MorpheusModification> GetNewVariableModificationPattern(int[] variableModificationArray, IEnumerable<KeyValuePair<int, UniqueModificationsList>> possibleVariableModifications)
        {
            var modification_pattern = new Dictionary<int, MorpheusModification>();

            foreach (KeyValuePair<int, UniqueModificationsList> kvp in possibleVariableModifications)
            {
                if (variableModificationArray[kvp.Key] > 0)
                {
                    modification_pattern.Add(kvp.Key, kvp.Value[variableModificationArray[kvp.Key] - 1]);
                }
            }

            return modification_pattern;
        }

        #endregion Private Methods
    }
}