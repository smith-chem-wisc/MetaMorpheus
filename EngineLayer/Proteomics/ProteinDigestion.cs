using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class ProteinDigestion
    {
        public Protease Protease { get; set; }
        public int MaximumMissedCleavages { get; set; }
        public IDigestionParams DigestionParams { get; set; }
        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; set; }
        public int? MinPeptidesLength { get; set; }
        public int? MaxPeptidesLength { get; set; }
        public IEnumerable<ModificationWithMass> AllKnownFixedModifications { get; set; }
        public List<ModificationWithMass> VariableModifications { get; set; }

        /// <summary>
        /// Initializes digestion object
        /// </summary>
        /// <param name="digestionParams"></param>
        /// <param name="allKnownFixedModifications"></param>
        /// <param name="variableModifications"></param>
        public ProteinDigestion(IDigestionParams digestionParams, IEnumerable<ModificationWithMass> allKnownFixedModifications, List<ModificationWithMass> variableModifications)
        {
            DigestionParams = digestionParams;
            Protease = digestionParams.Protease;
            MaximumMissedCleavages = digestionParams.MaxMissedCleavages;
            InitiatorMethionineBehavior = digestionParams.InitiatorMethionineBehavior;
            MinPeptidesLength = digestionParams.MinPeptideLength;
            MaxPeptidesLength = digestionParams.MaxPeptideLength;
            AllKnownFixedModifications = allKnownFixedModifications;
            VariableModifications = variableModifications;
        }

        /// <summary>
        /// Gets peptides for nonspecific digestion of a protein
        /// </summary>
        /// <param name="protein"></param>
        /// <returns></returns>
        public IEnumerable<PeptideWithSetModifications> NonSpecificDigestion(Protein protein)
        {
            List<ProteinInterval> intervals = new List<ProteinInterval>();
            List<int> oneBasedIndicesToCleaveAfter = Protease.GetDigestionSiteIndices(protein.BaseSequence);

            for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - MaximumMissedCleavages - 1; i++)
            {
                if (Protease.Retain(i, InitiatorMethionineBehavior, protein[0])
                    && Protease.OkayLength(oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1] - oneBasedIndicesToCleaveAfter[i], MinPeptidesLength, MaxPeptidesLength))
                {
                    intervals.Add(new ProteinInterval(oneBasedIndicesToCleaveAfter[i] + 1, oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1],
                        protein, oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1] - oneBasedIndicesToCleaveAfter[i], "semi"));
                }

                if (Protease.Cleave(i, InitiatorMethionineBehavior, protein[0])
                    && Protease.OkayLength(oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1] - 1, MinPeptidesLength, MaxPeptidesLength))
                {
                    intervals.Add(new ProteinInterval(2, oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1],
                        protein, oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1] - 1, "semi:M cleaved"));
                }
            }

            int lastIndex = oneBasedIndicesToCleaveAfter.Count - 1;
            int maxIndex = MaximumMissedCleavages < lastIndex ? MaximumMissedCleavages : lastIndex;
            for (int i = 1; i <= maxIndex; i++)
            {
                if (DigestionParams.TerminusTypeSemiProtease == TerminusType.N) //tricky, it's N because we want the extra peptide at the C terminus |_
                {
                    if (Protease.OkayLength(oneBasedIndicesToCleaveAfter[lastIndex] - oneBasedIndicesToCleaveAfter[lastIndex - i], MinPeptidesLength, MaxPeptidesLength))
                    {
                        intervals.Add(new ProteinInterval(oneBasedIndicesToCleaveAfter[lastIndex - i] + 1, oneBasedIndicesToCleaveAfter[lastIndex],
                            protein, oneBasedIndicesToCleaveAfter[lastIndex] - oneBasedIndicesToCleaveAfter[lastIndex - i], "semi"));
                    }
                }
                else //TerminusType.C
                {
                    if (Protease.OkayLength(oneBasedIndicesToCleaveAfter[i] - oneBasedIndicesToCleaveAfter[0], MinPeptidesLength, MaxPeptidesLength))
                    {
                        intervals.Add(new ProteinInterval(oneBasedIndicesToCleaveAfter[0] + 1, oneBasedIndicesToCleaveAfter[i],
                            protein, oneBasedIndicesToCleaveAfter[i] - oneBasedIndicesToCleaveAfter[0], "semi"));
                    }
                }
            }

            // Also digest using the proteolysis product start/end indices
            intervals.AddRange(
                protein.ProteolysisProducts
                    .Where(proteolysisProduct => proteolysisProduct.OneBasedBeginPosition != 1 || proteolysisProduct.OneBasedEndPosition != protein.Length)
                    .Select(proteolysisProduct => new ProteinInterval(proteolysisProduct.OneBasedBeginPosition.Value, proteolysisProduct.OneBasedEndPosition.Value,
                        protein, 0, proteolysisProduct.Type + " start")));

            return intervals.SelectMany(ivl => GetThePeptides(ivl, AllKnownFixedModifications, DigestionParams, VariableModifications));
        }

        /// <summary>
        /// Gets peptides for specific protease digestion of a protein
        /// </summary>
        /// <param name="protein"></param>
        /// <returns></returns>
        public IEnumerable<PeptideWithSetModifications> Digestion(Protein protein)
        {
            var intervals = Protease.GetDigestionIntervals(protein, MaximumMissedCleavages, InitiatorMethionineBehavior, MinPeptidesLength, MaxPeptidesLength);
            return intervals.SelectMany(ivl => GetThePeptides(ivl, AllKnownFixedModifications, DigestionParams, VariableModifications));
        }

        /// <summary>
        /// Gets the peptides for a specific protein interval
        /// </summary>
        /// <param name="interval"></param>
        /// <param name="allKnownFixedModifications"></param>
        /// <param name="digestionParams"></param>
        /// <param name="variableModifications"></param>
        /// <returns></returns>
        internal static IEnumerable<PeptideWithSetModifications> GetThePeptides(ProteinInterval interval,
            IEnumerable<ModificationWithMass> allKnownFixedModifications, IDigestionParams digestionParams, List<ModificationWithMass> variableModifications)
        {
            int peptideLength = interval.OneBasedEndResidueInProtein - interval.OneBasedStartResidueInProtein + 1;
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
                if (Gptmd.GptmdEngine.ModFits(variableModification, interval.Protein, 1, peptideLength, interval.OneBasedStartResidueInProtein)
                    && (variableModification.terminusLocalization == TerminusLocalization.NProt || variableModification.terminusLocalization == TerminusLocalization.NPep))
                {
                    pepNTermVariableMods.Add(variableModification);
                }

                for (int r = 0; r < peptideLength; r++)
                {
                    if (Gptmd.GptmdEngine.ModFits(variableModification, interval.Protein, r + 1, peptideLength, interval.OneBasedStartResidueInProtein + r)
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
                if (Gptmd.GptmdEngine.ModFits(variableModification, interval.Protein, peptideLength, peptideLength, interval.OneBasedStartResidueInProtein + peptideLength - 1)
                    && (variableModification.terminusLocalization == TerminusLocalization.ProtC || variableModification.terminusLocalization == TerminusLocalization.PepC))
                {
                    pepCTermVariableMods.Add(variableModification);
                }
            }

            // LOCALIZED MODS
            foreach (var kvp in interval.Protein.OneBasedPossibleLocalizedModifications)
            {
                bool inBounds = kvp.Key >= interval.OneBasedStartResidueInProtein && kvp.Key <= interval.OneBasedEndResidueInProtein;
                if (!inBounds)
                {
                    continue;
                }

                int locInPeptide = kvp.Key - interval.OneBasedStartResidueInProtein + 1;
                foreach (Modification modWithMass in kvp.Value)
                {
                    if (modWithMass is ModificationWithMass variableModification)
                    {
                        // Check if can be a n-term mod
                        if (locInPeptide == 1
                            && Gptmd.GptmdEngine.ModFits(variableModification, interval.Protein, 1, peptideLength, interval.OneBasedStartResidueInProtein)
                            && (variableModification.terminusLocalization == TerminusLocalization.NProt || variableModification.terminusLocalization == TerminusLocalization.NPep)
                            && !interval.Protein.IsDecoy)
                        {
                            pepNTermVariableMods.Add(variableModification);
                        }

                        for (int r = 0; r < peptideLength; r++)
                        {
                            if (locInPeptide == r + 1
                                && (interval.Protein.IsDecoy ||
                                (Gptmd.GptmdEngine.ModFits(variableModification, interval.Protein, r + 1, peptideLength, interval.OneBasedStartResidueInProtein + r)
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
                        if (locInPeptide == peptideLength
                            && Gptmd.GptmdEngine.ModFits(variableModification, interval.Protein, peptideLength, peptideLength, interval.OneBasedStartResidueInProtein + peptideLength - 1)
                            && (variableModification.terminusLocalization == TerminusLocalization.ProtC || variableModification.terminusLocalization == TerminusLocalization.PepC)
                            && !interval.Protein.IsDecoy)
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
                foreach (var ok in GetFixedModsOneIsNterminus(interval, peptideLength, allKnownFixedModifications))
                {
                    if (!kvp.ContainsKey(ok.Key))
                    {
                        numFixedMods++;
                        kvp.Add(ok.Key, ok.Value);
                    }
                }
                yield return new PeptideWithSetModifications(interval.Protein, interval.OneBasedStartResidueInProtein, interval.OneBasedEndResidueInProtein,
                    interval.PeptideString, interval.MissedCleavages, kvp, numFixedMods);
                variable_modification_isoforms++;
                if (variable_modification_isoforms == maximumVariableModificationIsoforms)
                {
                    yield break;
                }
            }
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

        private static Dictionary<int, ModificationWithMass> GetFixedModsOneIsNterminus(ProteinInterval interval, 
            int peptideLength, IEnumerable<ModificationWithMass> allKnownFixedModifications)
        {
            var fixedModsOneIsNterminus = new Dictionary<int, ModificationWithMass>(peptideLength + 3);
            foreach (ModificationWithMass mod in allKnownFixedModifications)
            {
                switch (mod.terminusLocalization)
                {
                    case TerminusLocalization.NProt:
                    case TerminusLocalization.NPep:
                        if (Gptmd.GptmdEngine.ModFits(mod, interval.Protein, 1, peptideLength, interval.OneBasedStartResidueInProtein))
                        {
                            fixedModsOneIsNterminus[1] = mod;
                        }
                        break;

                    case TerminusLocalization.Any:
                        for (int i = 2; i <= peptideLength + 1; i++)
                        {
                            if (Gptmd.GptmdEngine.ModFits(mod, interval.Protein, i - 1, peptideLength, interval.OneBasedStartResidueInProtein + i - 2))
                            {
                                fixedModsOneIsNterminus[i] = mod;
                            }
                        }
                        break;

                    case TerminusLocalization.PepC:
                    case TerminusLocalization.ProtC:
                        if (Gptmd.GptmdEngine.ModFits(mod, interval.Protein, peptideLength, peptideLength, interval.OneBasedStartResidueInProtein + peptideLength - 1))
                        {
                            fixedModsOneIsNterminus[peptideLength + 2] = mod;
                        }
                        break;
                }
            }
            return fixedModsOneIsNterminus;
        }
    }
}