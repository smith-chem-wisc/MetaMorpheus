using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class PeptideWithPossibleModifications : Peptide
    {

        #region Public Constructors

        public PeptideWithPossibleModifications(int oneBasedStartResidueNumberInProtein, int oneBasedEndResidueNumberInProtein, Protein parentProtein, int missedCleavages, string peptideDescription)
            : base(parentProtein, oneBasedStartResidueNumberInProtein, oneBasedEndResidueNumberInProtein)
        {
            OneBasedPossibleLocalizedModifications = Protein.OneBasedPossibleLocalizedModifications.Where(ok => ok.Key >= OneBasedStartResidueInProtein && ok.Key <= OneBasedEndResidueInProtein).ToDictionary(ok => ok.Key - OneBasedStartResidueInProtein + 1, ok => ok.Value);
            this.MissedCleavages = missedCleavages;
            this.PeptideDescription = peptideDescription;
        }

        #endregion Public Constructors

        #region Public Properties

        public Dictionary<int, List<MetaMorpheusModification>> OneBasedPossibleLocalizedModifications { get; private set; }
        public int MissedCleavages { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public IEnumerable<PeptideWithSetModifications> GetPeptideWithSetModifications(List<MetaMorpheusModification> variableModifications, int maximumVariableModificationIsoforms, int maxModsForPeptide, IEnumerable<MetaMorpheusModification> allKnownFixedModifications)
        {
            var two_based_possible_variable_and_localizeable_modifications = new Dictionary<int, UniqueModificationsCollection>(Length + 4);

            var pep_n_term_variable_mods = new UniqueModificationsCollection();
            two_based_possible_variable_and_localizeable_modifications.Add(1, pep_n_term_variable_mods);

            var pep_c_term_variable_mods = new UniqueModificationsCollection();
            two_based_possible_variable_and_localizeable_modifications.Add(Length + 2, pep_c_term_variable_mods);

            foreach (MetaMorpheusModification variable_modification in variableModifications)
            {
                if (variable_modification.ThisModificationType == ModificationType.ProteinNTerminus && (OneBasedStartResidueInProtein == 1 || (OneBasedStartResidueInProtein == 2 && Protein[0] == 'M'))
                    && (variable_modification.AminoAcid == char.MinValue || this[0] == variable_modification.AminoAcid))
                {
                    pep_n_term_variable_mods.Add(variable_modification);
                }

                if (variable_modification.ThisModificationType == ModificationType.PeptideNTerminus && (variable_modification.AminoAcid == char.MinValue || this[0] == variable_modification.AminoAcid))
                {
                    pep_n_term_variable_mods.Add(variable_modification);
                }

                for (int r = 0; r < Length; r++)
                {
                    if (variable_modification.ThisModificationType == ModificationType.AminoAcidResidue && this[r] == variable_modification.AminoAcid)
                    {
                        UniqueModificationsCollection residue_variable_mods;
                        if (!two_based_possible_variable_and_localizeable_modifications.TryGetValue(r + 2, out residue_variable_mods))
                        {
                            residue_variable_mods = new UniqueModificationsCollection();
                            residue_variable_mods.Add(variable_modification);
                            two_based_possible_variable_and_localizeable_modifications.Add(r + 2, residue_variable_mods);
                        }
                        else
                        {
                            residue_variable_mods.Add(variable_modification);
                        }
                    }
                }

                if (variable_modification.ThisModificationType == ModificationType.PeptideCTerminus && (variable_modification.AminoAcid == char.MinValue || this[Length - 1] == variable_modification.AminoAcid))
                {
                    pep_c_term_variable_mods.Add(variable_modification);
                }

                if (variable_modification.ThisModificationType == ModificationType.ProteinCTerminus && (OneBasedEndResidueInProtein == Protein.Length)
                    && (variable_modification.AminoAcid == char.MinValue || this[Length - 1] == variable_modification.AminoAcid))
                {
                    pep_c_term_variable_mods.Add(variable_modification);
                }
            }

            // LOCALIZED MODS
            if (OneBasedPossibleLocalizedModifications != null)
            {
                foreach (KeyValuePair<int, List<MetaMorpheusModification>> kvp in OneBasedPossibleLocalizedModifications)
                {
                    foreach (MetaMorpheusModification variable_modification in kvp.Value)
                    {
                        if (((variable_modification.ThisModificationType == ModificationType.ProteinNTerminus && !Protein.IsDecoy) ||
                            (variable_modification.ThisModificationType == ModificationType.ProteinCTerminus && Protein.IsDecoy))
                            && (OneBasedStartResidueInProtein == 1 || (OneBasedStartResidueInProtein == 2 && Protein[0] == 'M'))
                            && kvp.Key == 1)
                        {
                            pep_n_term_variable_mods.Add(variable_modification);
                        }

                        if (((variable_modification.ThisModificationType == ModificationType.PeptideNTerminus && !Protein.IsDecoy) ||
                            (variable_modification.ThisModificationType == ModificationType.PeptideCTerminus && Protein.IsDecoy))
                            && kvp.Key == 1)
                        {
                            pep_n_term_variable_mods.Add(variable_modification);
                        }

                        for (int r = 1; r <= Length; r++)
                        {
                            if (variable_modification.ThisModificationType == ModificationType.AminoAcidResidue
                                && kvp.Key == r)
                            {
                                UniqueModificationsCollection residue_variable_mods;
                                if (!two_based_possible_variable_and_localizeable_modifications.TryGetValue(r + 1, out residue_variable_mods))
                                {
                                    residue_variable_mods = new UniqueModificationsCollection();
                                    residue_variable_mods.Add(variable_modification);
                                    two_based_possible_variable_and_localizeable_modifications.Add(r + 1, residue_variable_mods);
                                }
                                else
                                {
                                    residue_variable_mods.Add(variable_modification);
                                }
                            }
                        }

                        if (((variable_modification.ThisModificationType == ModificationType.PeptideCTerminus && !Protein.IsDecoy) ||
                            (variable_modification.ThisModificationType == ModificationType.PeptideNTerminus && Protein.IsDecoy))
                            && kvp.Key == Length)
                        {
                            pep_c_term_variable_mods.Add(variable_modification);
                        }

                        if (((variable_modification.ThisModificationType == ModificationType.ProteinCTerminus && !Protein.IsDecoy) ||
                            (variable_modification.ThisModificationType == ModificationType.ProteinNTerminus && Protein.IsDecoy))
                            && (OneBasedEndResidueInProtein == Protein.Length)
                            && kvp.Key == Length)
                        {
                            pep_c_term_variable_mods.Add(variable_modification);
                        }
                    }
                }
            }

            int variable_modification_isoforms = 0;
            foreach (Dictionary<int, MetaMorpheusModification> kvp in GetVariableModificationPatterns(two_based_possible_variable_and_localizeable_modifications, maxModsForPeptide))
            {
                yield return new PeptideWithSetModifications(this, AddFixedMods(kvp, allKnownFixedModifications));
                variable_modification_isoforms++;
                if (variable_modification_isoforms == maximumVariableModificationIsoforms)
                    yield break;
            }
        }

        #endregion Public Methods

        #region Protected Methods

        protected IEnumerable<Dictionary<int, MetaMorpheusModification>> GetVariableModificationPatterns(Dictionary<int, UniqueModificationsCollection> possibleVariableModifications, int maxModsForPeptide)
        {
            if (possibleVariableModifications.Count == 0)
            {
                yield return null;
            }
            else
            {
                var possible_variable_modifications = new Dictionary<int, UniqueModificationsCollection>(possibleVariableModifications);

                int[] base_variable_modification_pattern = new int[Length + 4];
                var totalAvailableMods = possible_variable_modifications.Select(b => b.Value == null ? 0 : b.Value.Count).Sum();
                for (int variable_modifications = 0; variable_modifications <= Math.Min(totalAvailableMods, maxModsForPeptide); variable_modifications++)
                {
                    foreach (int[] variable_modification_pattern in GetVariableModificationPatterns(new List<KeyValuePair<int, UniqueModificationsCollection>>(possible_variable_modifications), possible_variable_modifications.Count - variable_modifications, base_variable_modification_pattern, 0))
                    {
                        yield return GetNewVariableModificationPattern(variable_modification_pattern, possible_variable_modifications);
                    }
                }
            }
        }

        #endregion Protected Methods

        #region Private Methods
        protected class UniqueModificationsCollection : List<MetaMorpheusModification>
        {

            #region Internal Methods

            internal new void Add(MetaMorpheusModification mod)
            {
                foreach (MetaMorpheusModification modHere in this)
                {
                    if (Math.Abs(modHere.PrecursorMassShift - mod.PrecursorMassShift) < 0.001)
                    {
                        return;
                    }
                }
                base.Add(mod);
            }

            #endregion Internal Methods
        }
        private static IEnumerable<int[]> GetVariableModificationPatterns(List<KeyValuePair<int, UniqueModificationsCollection>> possibleVariableModifications, int unmodifiedResiduesDesired, int[] variableModificationPattern, int index)
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

        private static Dictionary<int, MetaMorpheusModification> GetNewVariableModificationPattern(int[] variableModificationArray, IEnumerable<KeyValuePair<int, UniqueModificationsCollection>> possibleVariableModifications)
        {
            var modification_pattern = new Dictionary<int, MetaMorpheusModification>();

            foreach (KeyValuePair<int, UniqueModificationsCollection> kvp in possibleVariableModifications)
            {
                if (variableModificationArray[kvp.Key] > 0)
                {
                    modification_pattern.Add(kvp.Key, kvp.Value[variableModificationArray[kvp.Key] - 1]);
                }
            }

            return modification_pattern;
        }

        private Dictionary<int, MetaMorpheusModification> AddFixedMods(Dictionary<int, MetaMorpheusModification> allModsOneIsNterminus, IEnumerable<MetaMorpheusModification> allKnownFixedModifications)
        {
            MetaMorpheusModification val;
            for (int i = 0; i <= Length + 3; i++)
            {
                foreach (MetaMorpheusModification mod in allKnownFixedModifications)
                {
                    if (i == 0 && (OneBasedStartResidueInProtein == 1 || OneBasedStartResidueInProtein == 2))
                    {
                        if (mod.ThisModificationType == ModificationType.ProteinNTerminus && (mod.AminoAcid.Equals(this[0]) || mod.AminoAcid.Equals('\0')))
                        {
                            if (!allModsOneIsNterminus.TryGetValue(1, out val))
                                allModsOneIsNterminus.Add(1, mod);
                        }
                    }
                    else if (i == 1)
                    {
                        if (mod.ThisModificationType == ModificationType.PeptideNTerminus && (mod.AminoAcid.Equals(this[0]) || mod.AminoAcid.Equals('\0')))
                        {
                            if (!allModsOneIsNterminus.TryGetValue(1, out val))
                                allModsOneIsNterminus.Add(1, mod);
                        }
                    }
                    else if (i >= 2 && i <= Length + 1)
                    {
                        if (mod.ThisModificationType == ModificationType.AminoAcidResidue && (mod.AminoAcid.Equals(this[i - 2]) || mod.AminoAcid.Equals('\0')))
                        {
                            if (!allModsOneIsNterminus.TryGetValue(i, out val))
                                allModsOneIsNterminus.Add(i, mod);
                        }
                    }
                    else if (i == Length + 2)
                    {
                        if (mod.ThisModificationType == ModificationType.PeptideCTerminus && (mod.AminoAcid.Equals(this[Length - 1]) || mod.AminoAcid.Equals('\0')))
                        {
                            if (!allModsOneIsNterminus.TryGetValue(Length + 2, out val))
                                allModsOneIsNterminus.Add(Length + 2, mod);
                        }
                    }
                    else if (i == Length + 3 && OneBasedEndResidueInProtein == Protein.Length)
                    {
                        if (mod.ThisModificationType == ModificationType.ProteinCTerminus && (mod.AminoAcid.Equals(this[Length - 1]) || mod.AminoAcid.Equals('\0')))
                        {
                            if (!allModsOneIsNterminus.TryGetValue(Length + 2, out val))
                                allModsOneIsNterminus.Add(Length + 2, mod);
                        }
                    }
                }
            }
            return allModsOneIsNterminus;
        }

        #endregion Private Methods

    }
}