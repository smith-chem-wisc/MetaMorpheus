using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class PeptideWithPossibleModifications : Peptide
    {

        #region Private Fields

        private readonly Dictionary<int, ModificationWithMass> thisDictionaryOfFixedMods;

        #endregion Private Fields

        #region Internal Constructors

        internal PeptideWithPossibleModifications(int oneBasedStartResidueNumberInProtein, int oneBasedEndResidueNumberInProtein, Protein parentProtein, int missedCleavages, string peptideDescription, IEnumerable<ModificationWithMass> allKnownFixedModifications)
            : base(parentProtein, oneBasedStartResidueNumberInProtein, oneBasedEndResidueNumberInProtein)
        {
            this.MissedCleavages = missedCleavages;
            this.PeptideDescription = peptideDescription;
            thisDictionaryOfFixedMods = AddFixedMods(allKnownFixedModifications);
            this.numLocMods = this.Protein.OneBasedPossibleLocalizedModifications.Count(kvp => kvp.Key >= OneBasedStartResidueInProtein && kvp.Key <= OneBasedEndResidueInProtein);
        }

        #endregion Internal Constructors

        #region Public Properties

        public int MissedCleavages { get; private set; }
        public int numLocMods { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public IEnumerable<PeptideWithSetModifications> GetPeptidesWithSetModifications(List<ModificationWithMass> variableModifications, int maximumVariableModificationIsoforms, int maxModsForPeptide)
        {
            var two_based_possible_variable_and_localizeable_modifications = new Dictionary<int, UniqueModificationsCollection>(Length + 4);

            var pep_n_term_variable_mods = new UniqueModificationsCollection();
            two_based_possible_variable_and_localizeable_modifications.Add(1, pep_n_term_variable_mods);

            var pep_c_term_variable_mods = new UniqueModificationsCollection();
            two_based_possible_variable_and_localizeable_modifications.Add(Length + 2, pep_c_term_variable_mods);

            foreach (ModificationWithMass variable_modification in variableModifications)
            {
                // Check if can be a n-term mod
                if (Gptmd.GptmdEngine.ModFits(variable_modification, this.Protein, 1, this.Length, this.OneBasedStartResidueInProtein)
                    && (variable_modification.position == ModificationSites.NProt || variable_modification.position == ModificationSites.NPep))
                    pep_n_term_variable_mods.Add(variable_modification);

                for (int r = 0; r < Length; r++)
                {
                    if (Gptmd.GptmdEngine.ModFits(variable_modification, this.Protein, r + 1, this.Length, this.OneBasedStartResidueInProtein + r)
                        && variable_modification.position == ModificationSites.Any)
                    {
                        UniqueModificationsCollection residue_variable_mods;
                        if (!two_based_possible_variable_and_localizeable_modifications.TryGetValue(r + 2, out residue_variable_mods))
                        {
                            residue_variable_mods = new UniqueModificationsCollection();
                            residue_variable_mods.Add(variable_modification);
                            two_based_possible_variable_and_localizeable_modifications.Add(r + 2, residue_variable_mods);
                        }
                        else
                            residue_variable_mods.Add(variable_modification);
                    }
                }
                // Check if can be a c-term mod
                if (Gptmd.GptmdEngine.ModFits(variable_modification, this.Protein, Length, this.Length, this.OneBasedStartResidueInProtein + Length - 1)
                    && (variable_modification.position == ModificationSites.ProtC || variable_modification.position == ModificationSites.PepC))
                    pep_c_term_variable_mods.Add(variable_modification);
            }

            // LOCALIZED MODS
            foreach (var kvp in this.Protein.OneBasedPossibleLocalizedModifications)
            {
                if (kvp.Key >= OneBasedStartResidueInProtein && kvp.Key <= OneBasedEndResidueInProtein)
                {
                    int locInPeptide = kvp.Key - OneBasedStartResidueInProtein + 1;
                    foreach (Modification modMaybeWithMass in kvp.Value)
                    {
                        ModificationWithMass variable_modification = modMaybeWithMass as ModificationWithMass;
                        if (variable_modification != null)
                        {
                            // Check if can be a n-term mod
                            if (locInPeptide == 1
                                && Gptmd.GptmdEngine.ModFits(variable_modification, this.Protein, 1, this.Length, this.OneBasedStartResidueInProtein)
                                && (variable_modification.position == ModificationSites.NProt || variable_modification.position == ModificationSites.NPep))
                                pep_n_term_variable_mods.Add(variable_modification);

                            for (int r = 0; r < Length; r++)
                            {
                                if (locInPeptide == r + 1
                                    && Gptmd.GptmdEngine.ModFits(variable_modification, this.Protein, r + 1, this.Length, this.OneBasedStartResidueInProtein + r)
                                    && variable_modification.position == ModificationSites.Any)
                                {
                                    UniqueModificationsCollection residue_variable_mods;
                                    if (!two_based_possible_variable_and_localizeable_modifications.TryGetValue(r + 2, out residue_variable_mods))
                                    {
                                        residue_variable_mods = new UniqueModificationsCollection();
                                        residue_variable_mods.Add(variable_modification);
                                        two_based_possible_variable_and_localizeable_modifications.Add(r + 2, residue_variable_mods);
                                    }
                                    else
                                        residue_variable_mods.Add(variable_modification);
                                }
                            }
                            // Check if can be a c-term mod
                            if (locInPeptide == Length
                                && Gptmd.GptmdEngine.ModFits(variable_modification, this.Protein, Length, this.Length, this.OneBasedStartResidueInProtein + Length - 1)
                                && (variable_modification.position == ModificationSites.ProtC || variable_modification.position == ModificationSites.PepC))
                                pep_c_term_variable_mods.Add(variable_modification);
                        }
                    }
                }
            }

            int variable_modification_isoforms = 0;
            foreach (Dictionary<int, ModificationWithMass> kvp in GetVariableModificationPatterns(two_based_possible_variable_and_localizeable_modifications, maxModsForPeptide))
            {
                foreach (var ok in thisDictionaryOfFixedMods)
                    if (!kvp.ContainsKey(ok.Key))
                        kvp.Add(ok.Key, ok.Value);
                yield return new PeptideWithSetModifications(this, kvp);
                variable_modification_isoforms++;
                if (variable_modification_isoforms == maximumVariableModificationIsoforms)
                    yield break;
            }
        }

        #endregion Public Methods

        #region Protected Methods

        protected IEnumerable<Dictionary<int, ModificationWithMass>> GetVariableModificationPatterns(Dictionary<int, UniqueModificationsCollection> possibleVariableModifications, int maxModsForPeptide)
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

        private static Dictionary<int, ModificationWithMass> GetNewVariableModificationPattern(int[] variableModificationArray, IEnumerable<KeyValuePair<int, UniqueModificationsCollection>> possibleVariableModifications)
        {
            var modification_pattern = new Dictionary<int, ModificationWithMass>();

            foreach (KeyValuePair<int, UniqueModificationsCollection> kvp in possibleVariableModifications)
            {
                if (variableModificationArray[kvp.Key] > 0)
                {
                    modification_pattern.Add(kvp.Key, kvp.Value[variableModificationArray[kvp.Key] - 1]);
                }
            }

            return modification_pattern;
        }

        private Dictionary<int, ModificationWithMass> AddFixedMods(IEnumerable<ModificationWithMass> allKnownFixedModifications)
        {
            var fixedModsOneIsNterminus = new Dictionary<int, ModificationWithMass>(Length + 3);
            foreach (ModificationWithMass mod in allKnownFixedModifications)
            {
                switch (mod.position)
                {
                    case ModificationSites.NProt:
                    case ModificationSites.NPep:
                        if (Gptmd.GptmdEngine.ModFits(mod, this.Protein, 1, this.Length, this.OneBasedStartResidueInProtein))
                            fixedModsOneIsNterminus[1] = mod;
                        break;

                    case ModificationSites.Any:
                        for (int i = 2; i <= Length + 1; i++)
                            if (Gptmd.GptmdEngine.ModFits(mod, this.Protein, i - 1, this.Length, this.OneBasedStartResidueInProtein + i - 2))
                                fixedModsOneIsNterminus[i] = mod;
                        break;

                    case ModificationSites.PepC:
                    case ModificationSites.ProtC:
                        if (Gptmd.GptmdEngine.ModFits(mod, this.Protein, this.Length, this.Length, this.OneBasedStartResidueInProtein + this.Length - 1))
                            fixedModsOneIsNterminus[this.Length + 2] = mod;
                        break;
                }
            }
            return fixedModsOneIsNterminus;
        }

        #endregion Private Methods

        #region Protected Classes

        protected class UniqueModificationsCollection : List<ModificationWithMass>
        {

            #region Internal Methods

            internal new void Add(ModificationWithMass mod)
            {
                foreach (ModificationWithMass modHere in this)
                {
                    if (Math.Abs(modHere.monoisotopicMass - mod.monoisotopicMass) < 0.001 && Math.Abs(modHere.neutralLoss - mod.neutralLoss) < 0.001)
                    {
                        return;
                    }
                }
                base.Add(mod);
            }

            #endregion Internal Methods

        }

        #endregion Protected Classes

    }
}