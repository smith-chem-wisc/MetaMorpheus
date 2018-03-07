using Proteomics;
using System;
using System.Linq;
using System.Collections.Generic;

namespace EngineLayer
{
    public class Peptide
    {
        #region Private Fields

        private string baseSequence;

        #endregion Private Fields

        #region Protected Constructors

        internal Peptide(Protein protein, int oneBasedStartResidueInProtein, int oneBasedEndResidueInProtein, int? missedCleavages, string peptideDescription = null)
        {
            Protein = protein;
            OneBasedStartResidueInProtein = oneBasedStartResidueInProtein;
            OneBasedEndResidueInProtein = oneBasedEndResidueInProtein;
            Length = OneBasedEndResidueInProtein - OneBasedStartResidueInProtein + 1;
            MissedCleavages = missedCleavages;
            PeptideDescription = peptideDescription;
        }

        #endregion Protected Constructors

        #region Public Properties

        public Protein Protein { get; }
        public int OneBasedStartResidueInProtein { get; }
        public int OneBasedEndResidueInProtein { get; }
        public int? MissedCleavages { get; set; }
        public string PeptideDescription { get; }
        public int Length { get; }

        public virtual char PreviousAminoAcid
        {
            get
            {
                return OneBasedStartResidueInProtein > 1 ? Protein[OneBasedStartResidueInProtein - 2] : '-';
            }
        }

        public virtual char NextAminoAcid
        {
            get
            {
                return OneBasedEndResidueInProtein < Protein.Length ? Protein[OneBasedEndResidueInProtein] : '-';
            }
        }

        public string BaseSequence
        {
            get
            {
                if (baseSequence == null)
                    baseSequence = Protein.BaseSequence.Substring(OneBasedStartResidueInProtein - 1, Length);
                return baseSequence;
            }
        }

        #endregion Public Properties

        #region Public Indexers

        public char this[int zeroBasedIndex]
        {
            get
            {
                return Protein.BaseSequence[zeroBasedIndex + OneBasedStartResidueInProtein - 1];
            }
        }

        #endregion Public Indexers

        #region Modified Peptide Combinitorics

        /// <summary>
        /// Gets the peptides for a specific protein interval
        /// </summary>
        /// <param name="interval"></param>
        /// <param name="allKnownFixedModifications"></param>
        /// <param name="digestionParams"></param>
        /// <param name="variableModifications"></param>
        /// <returns></returns>
        internal IEnumerable<PeptideWithSetModifications> GetModifiedPeptides(IEnumerable<ModificationWithMass> allKnownFixedModifications, 
            IDigestionParams digestionParams, List<ModificationWithMass> variableModifications)
        {
            int peptideLength = OneBasedEndResidueInProtein - OneBasedStartResidueInProtein + 1;
            int maximumVariableModificationIsoforms = digestionParams.MaxModificationIsoforms;
            int maxModsForPeptide = digestionParams.MaxModsForPeptide;
            var twoBasedPossibleVariableAndLocalizeableModifications = new Dictionary<int, List<ModificationWithMass>>(peptideLength + 4);

            var pepNTermVariableMods = new List<ModificationWithMass>();
            twoBasedPossibleVariableAndLocalizeableModifications.Add(1, pepNTermVariableMods);

            var pepCTermVariableMods = new List<ModificationWithMass>();
            twoBasedPossibleVariableAndLocalizeableModifications.Add(peptideLength + 2, pepCTermVariableMods);

            foreach (ModificationWithMass variableModification in variableModifications)
            {
                // Check if can be a n-term mod
                if (CanBeNTerminalMod(variableModification, peptideLength))
                {
                    pepNTermVariableMods.Add(variableModification);
                }

                for (int r = 0; r < peptideLength; r++)
                {
                    if (Gptmd.GptmdEngine.ModFits(variableModification, Protein, r + 1, peptideLength, OneBasedStartResidueInProtein + r)
                        && variableModification.terminusLocalization == TerminusLocalization.Any)
                    {
                        if (!twoBasedPossibleVariableAndLocalizeableModifications.TryGetValue(r + 2, out List<ModificationWithMass> residueVariableMods))
                        {
                            residueVariableMods = new List<ModificationWithMass> { variableModification };
                            twoBasedPossibleVariableAndLocalizeableModifications.Add(r + 2, residueVariableMods);
                        }
                        else
                        {
                            residueVariableMods.Add(variableModification);
                        }
                    }
                }
                // Check if can be a c-term mod
                if (CanBeCTerminalMod(variableModification, peptideLength))
                {
                    pepCTermVariableMods.Add(variableModification);
                }
            }

            // LOCALIZED MODS
            foreach (var kvp in Protein.OneBasedPossibleLocalizedModifications)
            {
                bool inBounds = kvp.Key >= OneBasedStartResidueInProtein && kvp.Key <= OneBasedEndResidueInProtein;
                if (!inBounds)
                {
                    continue;
                }

                int locInPeptide = kvp.Key - OneBasedStartResidueInProtein + 1;
                foreach (Modification modWithMass in kvp.Value)
                {
                    if (modWithMass is ModificationWithMass variableModification)
                    {
                        // Check if can be a n-term mod
                        if (locInPeptide == 1 && CanBeNTerminalMod(variableModification, peptideLength) && !Protein.IsDecoy)
                        {
                            pepNTermVariableMods.Add(variableModification);
                        }

                        for (int r = 0; r < peptideLength; r++)
                        {
                            if (locInPeptide == r + 1
                                && (Protein.IsDecoy ||
                                (Gptmd.GptmdEngine.ModFits(variableModification, Protein, r + 1, peptideLength, OneBasedStartResidueInProtein + r)
                                    && variableModification.terminusLocalization == TerminusLocalization.Any)))
                            {
                                if (!twoBasedPossibleVariableAndLocalizeableModifications.TryGetValue(r + 2, out List<ModificationWithMass> residueVariableMods))
                                {
                                    residueVariableMods = new List<ModificationWithMass> { variableModification };
                                    twoBasedPossibleVariableAndLocalizeableModifications.Add(r + 2, residueVariableMods);
                                }
                                else
                                {
                                    residueVariableMods.Add(variableModification);
                                }
                            }
                        }
                        // Check if can be a c-term mod
                        if (locInPeptide == peptideLength && CanBeCTerminalMod(variableModification, peptideLength) && !Protein.IsDecoy)
                        {
                            pepCTermVariableMods.Add(variableModification);
                        }
                    }
                }
            }

            int variable_modification_isoforms = 0;

            foreach (Dictionary<int, ModificationWithMass> kvp in GetVariableModificationPatterns(twoBasedPossibleVariableAndLocalizeableModifications, maxModsForPeptide, peptideLength))
            {
                int numFixedMods = 0;
                foreach (var ok in GetFixedModsOneIsNterminus(peptideLength, allKnownFixedModifications))
                {
                    if (!kvp.ContainsKey(ok.Key))
                    {
                        numFixedMods++;
                        kvp.Add(ok.Key, ok.Value);
                    }
                }
                yield return new PeptideWithSetModifications(Protein, OneBasedStartResidueInProtein, OneBasedEndResidueInProtein,
                    PeptideDescription, MissedCleavages, kvp, numFixedMods);
                variable_modification_isoforms++;
                if (variable_modification_isoforms == maximumVariableModificationIsoforms)
                {
                    yield break;
                }
            }
        }

        /// <summary>
        /// Determines whether given modification can be an N-terminal modification
        /// </summary>
        /// <param name="variableModification"></param>
        /// <param name="peptideLength"></param>
        /// <returns></returns>
        private bool CanBeNTerminalMod(ModificationWithMass variableModification, int peptideLength)
        {
            return Gptmd.GptmdEngine.ModFits(variableModification, Protein, 1, peptideLength, OneBasedStartResidueInProtein)
                && (variableModification.terminusLocalization == TerminusLocalization.NProt || variableModification.terminusLocalization == TerminusLocalization.NPep);
        }

        /// <summary>
        /// Determines whether given modification can be a C-terminal modification
        /// </summary>
        /// <param name="variableModification"></param>
        /// <param name="peptideLength"></param>
        /// <returns></returns>
        private bool CanBeCTerminalMod(ModificationWithMass variableModification, int peptideLength)
        {
            return Gptmd.GptmdEngine.ModFits(variableModification, Protein, peptideLength, peptideLength, OneBasedStartResidueInProtein + peptideLength - 1)
                    && (variableModification.terminusLocalization == TerminusLocalization.ProtC || variableModification.terminusLocalization == TerminusLocalization.PepC);
        }

        private static IEnumerable<Dictionary<int, ModificationWithMass>> GetVariableModificationPatterns(Dictionary<int, List<ModificationWithMass>> possibleVariableModifications, int maxModsForPeptide, int peptideLength)
        {
            if (possibleVariableModifications.Count == 0)
            {
                yield return null;
            }
            else
            {
                var possible_variable_modifications = new Dictionary<int, List<ModificationWithMass>>(possibleVariableModifications);

                int[] base_variable_modification_pattern = new int[peptideLength + 4];
                var totalAvailableMods = possible_variable_modifications.Sum(b => b.Value == null ? 0 : b.Value.Count);
                for (int variable_modifications = 0; variable_modifications <= Math.Min(totalAvailableMods, maxModsForPeptide); variable_modifications++)
                {
                    foreach (int[] variable_modification_pattern in GetVariableModificationPatterns(new List<KeyValuePair<int, List<ModificationWithMass>>>(possible_variable_modifications),
                        possible_variable_modifications.Count - variable_modifications, base_variable_modification_pattern, 0))
                    {
                        yield return GetNewVariableModificationPattern(variable_modification_pattern, possible_variable_modifications);
                    }
                }
            }
        }

        private static IEnumerable<int[]> GetVariableModificationPatterns(List<KeyValuePair<int, List<ModificationWithMass>>> possibleVariableModifications,
            int unmodifiedResiduesDesired, int[] variableModificationPattern, int index)
        {
            if (index < possibleVariableModifications.Count - 1)
            {
                if (unmodifiedResiduesDesired > 0)
                {
                    variableModificationPattern[possibleVariableModifications[index].Key] = 0;
                    foreach (int[] new_variable_modification_pattern in GetVariableModificationPatterns(possibleVariableModifications,
                        unmodifiedResiduesDesired - 1, variableModificationPattern, index + 1))
                    {
                        yield return new_variable_modification_pattern;
                    }
                }
                if (unmodifiedResiduesDesired < possibleVariableModifications.Count - index)
                {
                    for (int i = 1; i <= possibleVariableModifications[index].Value.Count; i++)
                    {
                        variableModificationPattern[possibleVariableModifications[index].Key] = i;
                        foreach (int[] new_variable_modification_pattern in GetVariableModificationPatterns(possibleVariableModifications,
                            unmodifiedResiduesDesired, variableModificationPattern, index + 1))
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

        private static Dictionary<int, ModificationWithMass> GetNewVariableModificationPattern(int[] variableModificationArray,
            IEnumerable<KeyValuePair<int, List<ModificationWithMass>>> possibleVariableModifications)
        {
            var modification_pattern = new Dictionary<int, ModificationWithMass>();

            foreach (KeyValuePair<int, List<ModificationWithMass>> kvp in possibleVariableModifications)
            {
                if (variableModificationArray[kvp.Key] > 0)
                {
                    modification_pattern.Add(kvp.Key, kvp.Value[variableModificationArray[kvp.Key] - 1]);
                }
            }

            return modification_pattern;
        }

        private Dictionary<int, ModificationWithMass> GetFixedModsOneIsNterminus(int peptideLength, 
            IEnumerable<ModificationWithMass> allKnownFixedModifications)
        {
            var fixedModsOneIsNterminus = new Dictionary<int, ModificationWithMass>(peptideLength + 3);
            foreach (ModificationWithMass mod in allKnownFixedModifications)
            {
                switch (mod.terminusLocalization)
                {
                    case TerminusLocalization.NProt:
                    case TerminusLocalization.NPep:
                        if (Gptmd.GptmdEngine.ModFits(mod, Protein, 1, peptideLength, OneBasedStartResidueInProtein))
                        {
                            fixedModsOneIsNterminus[1] = mod;
                        }
                        break;

                    case TerminusLocalization.Any:
                        for (int i = 2; i <= peptideLength + 1; i++)
                        {
                            if (Gptmd.GptmdEngine.ModFits(mod, Protein, i - 1, peptideLength, OneBasedStartResidueInProtein + i - 2))
                            {
                                fixedModsOneIsNterminus[i] = mod;
                            }
                        }
                        break;

                    case TerminusLocalization.PepC:
                    case TerminusLocalization.ProtC:
                        if (Gptmd.GptmdEngine.ModFits(mod, Protein, peptideLength, peptideLength, OneBasedStartResidueInProtein + peptideLength - 1))
                        {
                            fixedModsOneIsNterminus[peptideLength + 2] = mod;
                        }
                        break;

                    default:
                        throw new NotSupportedException("This terminus localization is not supported.");
                }
            }
            return fixedModsOneIsNterminus;
        }

        #endregion Modified Peptide Combinitorics
    }
}