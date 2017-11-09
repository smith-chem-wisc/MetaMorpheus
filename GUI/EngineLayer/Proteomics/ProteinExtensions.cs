using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public static class ProteinExtensions
    {
        #region Public Methods

        public static IEnumerable<PeptideWithSetModifications> Digest(this Protein protein, DigestionParams digestionParams, IEnumerable<ModificationWithMass> allKnownFixedModifications, List<ModificationWithMass> variableModifications)
        {
            var protease = digestionParams.Protease;
            var maximumMissedCleavages = digestionParams.MaxMissedCleavages;
            var initiatorMethionineBehavior = digestionParams.InitiatorMethionineBehavior;
            var minPeptidesLength = digestionParams.MinPeptideLength;
            var maxPeptidesLength = digestionParams.MaxPeptideLength;
            List<int> oneBasedIndicesToCleaveAfter = new List<int>();

            if (digestionParams.SemiProteaseDigestion) //Used only for nonspecific search
            {  // these are the 1-based residue indices the protease cleaves AFTER
                oneBasedIndicesToCleaveAfter = protease.GetDigestionSiteIndices(protein.BaseSequence);
                // Cleave after 0, or before the first one
                oneBasedIndicesToCleaveAfter.Insert(0, 0);
                // Cleave after Length index
                oneBasedIndicesToCleaveAfter.Add(protein.Length);
                for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - maximumMissedCleavages - 1; i++)
                {
                    // Retain!
                    if ((i != 0 || initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave || protein[0] != 'M')
                        && (!minPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1] - oneBasedIndicesToCleaveAfter[i] >= minPeptidesLength)
                        && (!maxPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1] - oneBasedIndicesToCleaveAfter[i] <= maxPeptidesLength))
                        foreach (PeptideWithSetModifications p in GetThePeptides(oneBasedIndicesToCleaveAfter[i] + 1, oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1], protein, oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1] - oneBasedIndicesToCleaveAfter[i], "semi", allKnownFixedModifications, digestionParams, variableModifications))
                            yield return p;

                    // Cleave!
                    if ((i == 0 && initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && protein[0] == 'M')
                        && (!minPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1] - 1 >= minPeptidesLength)
                        && (!maxPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1] - 1 <= maxPeptidesLength))
                        foreach (PeptideWithSetModifications p in GetThePeptides(2, oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1], protein, oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1] - 1, "semi:M cleaved", allKnownFixedModifications, digestionParams, variableModifications))
                            yield return p;
                }

                int lastIndex = oneBasedIndicesToCleaveAfter.Count - 1;
                int maxIndex = maximumMissedCleavages < lastIndex ? maximumMissedCleavages : lastIndex;
                if (digestionParams.TerminusTypeSemiProtease == TerminusType.N) //tricky, it's N because we want the extra peptide at the C terminus |_
                {
                    for (int i = 1; i <= maxIndex; i++)
                    {
                        if ((!minPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[lastIndex] - oneBasedIndicesToCleaveAfter[lastIndex - i] >= minPeptidesLength) &&
                            (!maxPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[lastIndex] - oneBasedIndicesToCleaveAfter[lastIndex - i] <= maxPeptidesLength))
                        {
                            foreach (PeptideWithSetModifications p in GetThePeptides(oneBasedIndicesToCleaveAfter[lastIndex - i] + 1, oneBasedIndicesToCleaveAfter[lastIndex], protein, oneBasedIndicesToCleaveAfter[lastIndex] - oneBasedIndicesToCleaveAfter[lastIndex - i], "semi", allKnownFixedModifications, digestionParams, variableModifications))
                                yield return p;
                        }
                    }
                }
                else //TerminusType.C
                {
                    for (int i = 1; i <= maxIndex; i++)
                    {
                        if ((!minPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i] - oneBasedIndicesToCleaveAfter[0] >= minPeptidesLength) &&
                            (!maxPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i] - oneBasedIndicesToCleaveAfter[0] <= maxPeptidesLength))
                        {
                            foreach (PeptideWithSetModifications p in GetThePeptides(oneBasedIndicesToCleaveAfter[0] + 1, oneBasedIndicesToCleaveAfter[i], protein, oneBasedIndicesToCleaveAfter[i] - oneBasedIndicesToCleaveAfter[0], "semi", allKnownFixedModifications, digestionParams, variableModifications))
                                yield return p;
                        }
                    }
                }
                // Also digest using the proteolysis product start/end indices
                foreach (var proteolysisProduct in protein.ProteolysisProducts)
                    if (proteolysisProduct.OneBasedBeginPosition != 1 || proteolysisProduct.OneBasedEndPosition != protein.Length)
                    {
                        foreach (PeptideWithSetModifications p in GetThePeptides(proteolysisProduct.OneBasedBeginPosition.Value, proteolysisProduct.OneBasedEndPosition.Value, protein, 0, proteolysisProduct.Type + " start", allKnownFixedModifications, digestionParams, variableModifications))
                            yield return p;
                    }
            }
            else
            {
                switch (protease.CleavageSpecificity)
                {
                    case CleavageSpecificity.Full:
                        // these are the 1-based residue indices the protease cleaves AFTER
                        oneBasedIndicesToCleaveAfter = protease.GetDigestionSiteIndices(protein.BaseSequence);
                        // Cleave after 0, or before the first one
                        oneBasedIndicesToCleaveAfter.Insert(0, 0);
                        // Cleave after Length index
                        oneBasedIndicesToCleaveAfter.Add(protein.Length);
                        for (int missed_cleavages = 0; missed_cleavages <= maximumMissedCleavages; missed_cleavages++)
                        {
                            for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - missed_cleavages - 1; i++)
                            {
                                // Retain!
                                if (i != 0 || initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave || protein[0] != 'M')
                                {
                                    if ((!minPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1] - oneBasedIndicesToCleaveAfter[i] >= minPeptidesLength) &&
                                        (!maxPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1] - oneBasedIndicesToCleaveAfter[i] <= maxPeptidesLength))
                                    {
                                        foreach (PeptideWithSetModifications p in GetThePeptides(oneBasedIndicesToCleaveAfter[i] + 1, oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1], protein, missed_cleavages, "full", allKnownFixedModifications, digestionParams, variableModifications))
                                            yield return p;
                                    }
                                }
                                // Cleave!
                                if (i == 0 && initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && protein[0] == 'M')
                                {
                                    if ((!minPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1] - 1 >= minPeptidesLength) &&
                                        (!maxPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1] - 1 <= maxPeptidesLength))
                                    {
                                        foreach (PeptideWithSetModifications p in GetThePeptides(2, oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1], protein, missed_cleavages, "full:M cleaved", allKnownFixedModifications, digestionParams, variableModifications))
                                            yield return p;
                                    }
                                }
                            }

                            // Also digest using the proteolysis product start/end indices
                            foreach (var proteolysisProduct in protein.ProteolysisProducts)
                                if (proteolysisProduct.OneBasedBeginPosition != 1 || proteolysisProduct.OneBasedEndPosition != protein.Length)
                                {
                                    int i = 0;
                                    while (oneBasedIndicesToCleaveAfter[i] < proteolysisProduct.OneBasedBeginPosition)
                                        i++;
                                    // Start peptide
                                    if (i + missed_cleavages < oneBasedIndicesToCleaveAfter.Count && oneBasedIndicesToCleaveAfter[i + missed_cleavages] <= proteolysisProduct.OneBasedEndPosition && proteolysisProduct.OneBasedBeginPosition.HasValue)
                                    {
                                        if ((!minPeptidesLength.HasValue || (oneBasedIndicesToCleaveAfter[i + missed_cleavages] - proteolysisProduct.OneBasedBeginPosition.Value + 1) >= minPeptidesLength) &&
                                            (!maxPeptidesLength.HasValue || (oneBasedIndicesToCleaveAfter[i + missed_cleavages] - proteolysisProduct.OneBasedBeginPosition.Value + 1) <= maxPeptidesLength))
                                        {
                                            foreach (PeptideWithSetModifications p in GetThePeptides(proteolysisProduct.OneBasedBeginPosition.Value, oneBasedIndicesToCleaveAfter[i + missed_cleavages], protein, missed_cleavages, proteolysisProduct.Type + " start", allKnownFixedModifications, digestionParams, variableModifications))
                                                yield return p;
                                        }
                                    }
                                    while (oneBasedIndicesToCleaveAfter[i] < proteolysisProduct.OneBasedEndPosition)
                                        i++;
                                    // End
                                    if (i - missed_cleavages - 1 >= 0 && oneBasedIndicesToCleaveAfter[i - missed_cleavages - 1] + 1 >= proteolysisProduct.OneBasedBeginPosition && proteolysisProduct.OneBasedEndPosition.HasValue)
                                    {
                                        if ((!minPeptidesLength.HasValue || (proteolysisProduct.OneBasedEndPosition.Value - oneBasedIndicesToCleaveAfter[i - missed_cleavages - 1]) >= minPeptidesLength) &&
                                            (!maxPeptidesLength.HasValue || (proteolysisProduct.OneBasedEndPosition.Value - oneBasedIndicesToCleaveAfter[i - missed_cleavages - 1]) <= maxPeptidesLength))
                                        {
                                            foreach (PeptideWithSetModifications p in GetThePeptides(oneBasedIndicesToCleaveAfter[i - missed_cleavages - 1] + 1, proteolysisProduct.OneBasedEndPosition.Value, protein, missed_cleavages, proteolysisProduct.Type + " end", allKnownFixedModifications, digestionParams, variableModifications))
                                                yield return p;
                                        }
                                    }
                                }
                        }
                        break;

                    case CleavageSpecificity.SingleN:
                        //cleave in one spot
                        for (int index = 1; index <= protein.Length; index++) //position BEFORE the amino acid
                        {
                            if ((!minPeptidesLength.HasValue || (protein.Length - index + 1) >= minPeptidesLength)) //&&
                            {
                                foreach (PeptideWithSetModifications p in GetThePeptides(index, Math.Min(protein.Length, index + (maxPeptidesLength ?? 50)), protein, 0, "SingleN", allKnownFixedModifications, digestionParams, variableModifications))
                                    yield return p;
                            }
                        }
                        break;

                    case CleavageSpecificity.SingleC:
                        //cleave in one spot
                        for (int index = 1; index <= protein.Length; index++) //position BEFORE the amino acid
                        {
                            if ((!minPeptidesLength.HasValue || (index) >= minPeptidesLength))//&&
                            {
                                foreach (PeptideWithSetModifications p in GetThePeptides(Math.Max(1, index - (maxPeptidesLength ?? 50)), index, protein, 0, "SingleC", allKnownFixedModifications, digestionParams, variableModifications))
                                    yield return p;
                            }
                        }
                        break;

                    case CleavageSpecificity.Semi:
                        // these are the 1-based residue indices the protease cleaves AFTER
                        oneBasedIndicesToCleaveAfter = protease.GetDigestionSiteIndices(protein.BaseSequence);
                        // Cleave after 0, or before the first one
                        oneBasedIndicesToCleaveAfter.Insert(0, 0);
                        // Cleave after Length index
                        oneBasedIndicesToCleaveAfter.Add(protein.Length);
                        //Main Body
                        for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - maximumMissedCleavages - 1; i++)
                        {
                            // Retain!
                            if (i != 0 || initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave || protein[0] != 'M')
                            {
                                int nTermIndexofProtein = oneBasedIndicesToCleaveAfter[i];
                                int cTermIndexofProtein = oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1];

                                if ((!minPeptidesLength.HasValue || cTermIndexofProtein - nTermIndexofProtein >= minPeptidesLength) &&
                                    (!maxPeptidesLength.HasValue || cTermIndexofProtein - nTermIndexofProtein <= maxPeptidesLength))
                                {
                                    foreach (PeptideWithSetModifications p in GetThePeptides(nTermIndexofProtein + 1, cTermIndexofProtein, protein, cTermIndexofProtein - nTermIndexofProtein, "semi", allKnownFixedModifications, digestionParams, variableModifications))
                                        yield return p;
                                }
                                //fixedC
                                for (int j = nTermIndexofProtein + 1; j < cTermIndexofProtein; j++)
                                {
                                    if ((!minPeptidesLength.HasValue || cTermIndexofProtein - j >= minPeptidesLength) &&
                                        (!maxPeptidesLength.HasValue || cTermIndexofProtein - j <= maxPeptidesLength))
                                    {
                                        foreach (PeptideWithSetModifications p in GetThePeptides(j + 1, cTermIndexofProtein, protein, cTermIndexofProtein - j, "semiC", allKnownFixedModifications, digestionParams, variableModifications))
                                            yield return p;
                                    }
                                }
                                //fixedN
                                for (int j = cTermIndexofProtein - 1; j > nTermIndexofProtein; j--)
                                {
                                    if ((!minPeptidesLength.HasValue || j - nTermIndexofProtein >= minPeptidesLength) &&
                                        (!maxPeptidesLength.HasValue || j - nTermIndexofProtein <= maxPeptidesLength))
                                    {
                                        foreach (PeptideWithSetModifications p in GetThePeptides(nTermIndexofProtein + 1, j, protein, j - nTermIndexofProtein, "semiN", allKnownFixedModifications, digestionParams, variableModifications))
                                            yield return p;
                                    }
                                }
                            }
                            // Cleave!
                            if (i == 0 && initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && protein[0] == 'M')
                            {
                                int nTermIndexOfProtein = 1;
                                int cTermIndexOfProtein = oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1];
                                if ((!minPeptidesLength.HasValue || cTermIndexOfProtein - nTermIndexOfProtein >= minPeptidesLength) &&
                                    (!maxPeptidesLength.HasValue || cTermIndexOfProtein - nTermIndexOfProtein <= maxPeptidesLength))
                                {
                                    foreach (PeptideWithSetModifications p in GetThePeptides(nTermIndexOfProtein + 1, cTermIndexOfProtein, protein, cTermIndexOfProtein - nTermIndexOfProtein, "semi:M cleaved", allKnownFixedModifications, digestionParams, variableModifications))
                                        yield return p;
                                }
                                //fixedC
                                for (int j = nTermIndexOfProtein + 1; j < cTermIndexOfProtein; j++)
                                {
                                    if ((!minPeptidesLength.HasValue || cTermIndexOfProtein - j >= minPeptidesLength) &&
                                        (!maxPeptidesLength.HasValue || cTermIndexOfProtein - j <= maxPeptidesLength))
                                    {
                                        foreach (PeptideWithSetModifications p in GetThePeptides(j + 1, cTermIndexOfProtein, protein, cTermIndexOfProtein - j, "semiC", allKnownFixedModifications, digestionParams, variableModifications))
                                            yield return p;
                                    }
                                }
                                //fixedN
                                for (int j = cTermIndexOfProtein - 1; j > nTermIndexOfProtein; j--)
                                {
                                    if ((!minPeptidesLength.HasValue || j - nTermIndexOfProtein >= minPeptidesLength) &&
                                        (!maxPeptidesLength.HasValue || j - nTermIndexOfProtein <= maxPeptidesLength))
                                    {
                                        foreach (PeptideWithSetModifications p in GetThePeptides(nTermIndexOfProtein + 1, j, protein, j - nTermIndexOfProtein, "semiN", allKnownFixedModifications, digestionParams, variableModifications))
                                            yield return p;
                                    }
                                }
                            }
                        }
                        int lastIndexSemi = oneBasedIndicesToCleaveAfter.Count - 1;
                        int maxIndexSemi = maximumMissedCleavages < lastIndexSemi ? maximumMissedCleavages : lastIndexSemi;
                        //Fringe C-term peptides
                        for (int i = 1; i <= maxIndexSemi; i++)
                        {
                            //fixedN
                            int nTermIndexOfProtein = oneBasedIndicesToCleaveAfter[lastIndexSemi - i];
                            int cTermIndexOfProtein = oneBasedIndicesToCleaveAfter[lastIndexSemi] - 1;
                            for (int j = cTermIndexOfProtein; j > nTermIndexOfProtein; j--)
                            {
                                if ((!minPeptidesLength.HasValue || j - nTermIndexOfProtein >= minPeptidesLength) &&
                                    (!maxPeptidesLength.HasValue || j - nTermIndexOfProtein <= maxPeptidesLength))
                                {
                                    foreach (PeptideWithSetModifications p in GetThePeptides(nTermIndexOfProtein + 1, j, protein, j - nTermIndexOfProtein, "semiN", allKnownFixedModifications, digestionParams, variableModifications))
                                        yield return p;
                                }
                            }
                        }
                        //Fringe N-term peptides
                        for (int i = 1; i <= maxIndexSemi; i++)
                        {
                            //fixedC
                            int nTermIndexOfProtein = oneBasedIndicesToCleaveAfter[0];
                            int cTermIndexOfProtein = oneBasedIndicesToCleaveAfter[i];
                            for (int j = nTermIndexOfProtein + 1; j < cTermIndexOfProtein; j++)
                            {
                                if ((!minPeptidesLength.HasValue || cTermIndexOfProtein - j >= minPeptidesLength) &&
                                    (!maxPeptidesLength.HasValue || cTermIndexOfProtein - j <= maxPeptidesLength))
                                {
                                    foreach (PeptideWithSetModifications p in GetThePeptides(j + 1, cTermIndexOfProtein, protein, cTermIndexOfProtein - j, "semiC", allKnownFixedModifications, digestionParams, variableModifications))
                                        yield return p;
                                }
                            }
                        }

                        // Also digest using the proteolysis product start/end indices
                        foreach (var proteolysisProduct in protein.ProteolysisProducts)
                            if (proteolysisProduct.OneBasedBeginPosition != 1 || proteolysisProduct.OneBasedEndPosition != protein.Length)
                            {
                                int i = 0;
                                while (oneBasedIndicesToCleaveAfter[i] < proteolysisProduct.OneBasedBeginPosition)//"<" to prevent additions if same index as residues
                                    i++;
                                // Start peptide
                                for (int j = proteolysisProduct.OneBasedBeginPosition.Value; j < oneBasedIndicesToCleaveAfter[i]; j++)
                                {
                                    if ((!minPeptidesLength.HasValue || j - proteolysisProduct.OneBasedBeginPosition + 1 >= minPeptidesLength) &&
                                        (!maxPeptidesLength.HasValue || j - proteolysisProduct.OneBasedBeginPosition + 1 <= maxPeptidesLength))
                                    {
                                        foreach (PeptideWithSetModifications p in GetThePeptides(proteolysisProduct.OneBasedBeginPosition.Value, j, protein, j - proteolysisProduct.OneBasedBeginPosition.Value, proteolysisProduct.Type + " start", allKnownFixedModifications, digestionParams, variableModifications))
                                            yield return p;
                                    }
                                }
                                while (oneBasedIndicesToCleaveAfter[i] <= proteolysisProduct.OneBasedEndPosition) //"<=" to prevent additions if same index as residues, since i-- is below
                                    i++;
                                i--;
                                // End
                                for (int j = oneBasedIndicesToCleaveAfter[i] + 1; j < proteolysisProduct.OneBasedEndPosition.Value; j++)
                                {
                                    if ((!minPeptidesLength.HasValue || proteolysisProduct.OneBasedEndPosition - j + 1 >= minPeptidesLength) &&
                                        (!maxPeptidesLength.HasValue || proteolysisProduct.OneBasedEndPosition - j + 1 <= maxPeptidesLength))
                                    {
                                        foreach (PeptideWithSetModifications p in GetThePeptides(j, proteolysisProduct.OneBasedEndPosition.Value, protein, proteolysisProduct.OneBasedEndPosition.Value - j, proteolysisProduct.Type + " end", allKnownFixedModifications, digestionParams, variableModifications))
                                            yield return p;
                                    }
                                }
                            }
                        break;

                    case CleavageSpecificity.None:
                        if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave || protein[0] != 'M')
                        {
                            if ((!minPeptidesLength.HasValue || protein.Length >= minPeptidesLength) &&
                                            (!maxPeptidesLength.HasValue || protein.Length <= maxPeptidesLength))
                            {
                                foreach (PeptideWithSetModifications p in GetThePeptides(1, protein.Length, protein, 0, "full", allKnownFixedModifications, digestionParams, variableModifications))
                                    yield return p;
                            }
                        }
                        if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && protein[0] == 'M')
                        {
                            if ((!minPeptidesLength.HasValue || protein.Length - 1 >= minPeptidesLength) &&
                                            (!maxPeptidesLength.HasValue || protein.Length - 1 <= maxPeptidesLength))
                            {
                                foreach (PeptideWithSetModifications p in GetThePeptides(2, protein.Length, protein, 0, "full:M cleaved", allKnownFixedModifications, digestionParams, variableModifications))
                                    yield return p;
                            }
                        }

                        // Also digest using the proteolysis product start/end indices
                        foreach (var proteolysisProduct in protein.ProteolysisProducts)
                            if (proteolysisProduct.OneBasedEndPosition.HasValue && proteolysisProduct.OneBasedBeginPosition.HasValue)
                            {
                                if ((!minPeptidesLength.HasValue || (proteolysisProduct.OneBasedEndPosition.Value - proteolysisProduct.OneBasedBeginPosition.Value + 1) >= minPeptidesLength) &&
                                               (!maxPeptidesLength.HasValue || (proteolysisProduct.OneBasedEndPosition.Value - proteolysisProduct.OneBasedBeginPosition.Value + 1) <= maxPeptidesLength))
                                {
                                    foreach (PeptideWithSetModifications p in GetThePeptides(proteolysisProduct.OneBasedBeginPosition.Value, proteolysisProduct.OneBasedEndPosition.Value, protein, 0, proteolysisProduct.Type, allKnownFixedModifications, digestionParams, variableModifications))
                                        yield return p;
                                }
                            }
                        break;

                    default:
                        throw new NotImplementedException();
                }
            }
        }

        #endregion Public Methods

        #region Private Methods

        private static IEnumerable<PeptideWithSetModifications> GetThePeptides(int OneBasedStartResidueInProtein, int OneBasedEndResidueInProtein, Protein protein, int missedCleavages, string peptideString, IEnumerable<ModificationWithMass> allKnownFixedModifications, DigestionParams digestionParams, List<ModificationWithMass> variableModifications)
        {
            int peptideLength = OneBasedEndResidueInProtein - OneBasedStartResidueInProtein + 1;
            int maximumVariableModificationIsoforms = digestionParams.MaxModificationIsoforms;
            int maxModsForPeptide = digestionParams.MaxModsForPeptide;
            var two_based_possible_variable_and_localizeable_modifications = new Dictionary<int, List<ModificationWithMass>>(peptideLength + 4);

            var pep_n_term_variable_mods = new List<ModificationWithMass>();
            two_based_possible_variable_and_localizeable_modifications.Add(1, pep_n_term_variable_mods);

            var pep_c_term_variable_mods = new List<ModificationWithMass>();
            two_based_possible_variable_and_localizeable_modifications.Add(peptideLength + 2, pep_c_term_variable_mods);

            foreach (ModificationWithMass variable_modification in variableModifications)
            {
                // Check if can be a n-term mod
                if (Gptmd.GptmdEngine.ModFits(variable_modification, protein, 1, peptideLength, OneBasedStartResidueInProtein)
                    && (variable_modification.terminusLocalization == TerminusLocalization.NProt || variable_modification.terminusLocalization == TerminusLocalization.NPep))
                    pep_n_term_variable_mods.Add(variable_modification);

                for (int r = 0; r < peptideLength; r++)
                {
                    if (Gptmd.GptmdEngine.ModFits(variable_modification, protein, r + 1, peptideLength, OneBasedStartResidueInProtein + r)
                        && variable_modification.terminusLocalization == TerminusLocalization.Any)
                    {
                        if (!two_based_possible_variable_and_localizeable_modifications.TryGetValue(r + 2, out List<ModificationWithMass> residue_variable_mods))
                        {
                            residue_variable_mods = new List<ModificationWithMass>
                            {
                                variable_modification
                            };
                            two_based_possible_variable_and_localizeable_modifications.Add(r + 2, residue_variable_mods);
                        }
                        else
                            residue_variable_mods.Add(variable_modification);
                    }
                }
                // Check if can be a c-term mod
                if (Gptmd.GptmdEngine.ModFits(variable_modification, protein, peptideLength, peptideLength, OneBasedStartResidueInProtein + peptideLength - 1)
                    && (variable_modification.terminusLocalization == TerminusLocalization.ProtC || variable_modification.terminusLocalization == TerminusLocalization.PepC))
                    pep_c_term_variable_mods.Add(variable_modification);
            }

            // LOCALIZED MODS
            foreach (var kvp in protein.OneBasedPossibleLocalizedModifications)
            {
                if (kvp.Key >= OneBasedStartResidueInProtein && kvp.Key <= OneBasedEndResidueInProtein)
                {
                    int locInPeptide = kvp.Key - OneBasedStartResidueInProtein + 1;
                    foreach (Modification modMaybeWithMass in kvp.Value)
                    {
                        if (modMaybeWithMass is ModificationWithMass variable_modification)
                        {
                            // Check if can be a n-term mod
                            if (locInPeptide == 1
                                && Gptmd.GptmdEngine.ModFits(variable_modification, protein, 1, peptideLength, OneBasedStartResidueInProtein)
                                && (variable_modification.terminusLocalization == TerminusLocalization.NProt || variable_modification.terminusLocalization == TerminusLocalization.NPep)
                                && !protein.IsDecoy)
                                pep_n_term_variable_mods.Add(variable_modification);

                            for (int r = 0; r < peptideLength; r++)
                            {
                                if (locInPeptide == r + 1
                                    && (protein.IsDecoy || (Gptmd.GptmdEngine.ModFits(variable_modification, protein, r + 1, peptideLength, OneBasedStartResidueInProtein + r)
                                    && variable_modification.terminusLocalization == TerminusLocalization.Any)))
                                {
                                    if (!two_based_possible_variable_and_localizeable_modifications.TryGetValue(r + 2, out List<ModificationWithMass> residue_variable_mods))
                                    {
                                        residue_variable_mods = new List<ModificationWithMass>
                                        {
                                            variable_modification
                                        };
                                        two_based_possible_variable_and_localizeable_modifications.Add(r + 2, residue_variable_mods);
                                    }
                                    else
                                        residue_variable_mods.Add(variable_modification);
                                }
                            }
                            // Check if can be a c-term mod
                            if (locInPeptide == peptideLength
                                && Gptmd.GptmdEngine.ModFits(variable_modification, protein, peptideLength, peptideLength, OneBasedStartResidueInProtein + peptideLength - 1)
                                && (variable_modification.terminusLocalization == TerminusLocalization.ProtC || variable_modification.terminusLocalization == TerminusLocalization.PepC)
                                && !protein.IsDecoy)
                                pep_c_term_variable_mods.Add(variable_modification);
                        }
                    }
                }
            }

            int variable_modification_isoforms = 0;

            var fixedModsOneIsNterminus = new Dictionary<int, ModificationWithMass>(peptideLength + 3);
            foreach (ModificationWithMass mod in allKnownFixedModifications)
            {
                switch (mod.terminusLocalization)
                {
                    case TerminusLocalization.NProt:
                    case TerminusLocalization.NPep:
                        if (Gptmd.GptmdEngine.ModFits(mod, protein, 1, peptideLength, OneBasedStartResidueInProtein))
                            fixedModsOneIsNterminus[1] = mod;
                        break;

                    case TerminusLocalization.Any:
                        for (int i = 2; i <= peptideLength + 1; i++)
                            if (Gptmd.GptmdEngine.ModFits(mod, protein, i - 1, peptideLength, OneBasedStartResidueInProtein + i - 2))
                                fixedModsOneIsNterminus[i] = mod;
                        break;

                    case TerminusLocalization.PepC:
                    case TerminusLocalization.ProtC:
                        if (Gptmd.GptmdEngine.ModFits(mod, protein, peptideLength, peptideLength, OneBasedStartResidueInProtein + peptideLength - 1))
                            fixedModsOneIsNterminus[peptideLength + 2] = mod;
                        break;
                }
            }

            foreach (Dictionary<int, ModificationWithMass> kvp in GetVariableModificationPatterns(two_based_possible_variable_and_localizeable_modifications, maxModsForPeptide, peptideLength))
            {
                int numFixedMods = 0;
                foreach (var ok in fixedModsOneIsNterminus)
                    if (!kvp.ContainsKey(ok.Key))
                    {
                        numFixedMods++;
                        kvp.Add(ok.Key, ok.Value);
                    }
                yield return new PeptideWithSetModifications(protein, OneBasedStartResidueInProtein, OneBasedEndResidueInProtein, peptideString, missedCleavages, kvp, numFixedMods);
                variable_modification_isoforms++;
                if (variable_modification_isoforms == maximumVariableModificationIsoforms)
                    yield break;
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
                var totalAvailableMods = possible_variable_modifications.Select(b => b.Value == null ? 0 : b.Value.Count).Sum();
                for (int variable_modifications = 0; variable_modifications <= Math.Min(totalAvailableMods, maxModsForPeptide); variable_modifications++)
                {
                    foreach (int[] variable_modification_pattern in GetVariableModificationPatterns(new List<KeyValuePair<int, List<ModificationWithMass>>>(possible_variable_modifications), possible_variable_modifications.Count - variable_modifications, base_variable_modification_pattern, 0))
                    {
                        yield return GetNewVariableModificationPattern(variable_modification_pattern, possible_variable_modifications);
                    }
                }
            }
        }

        private static IEnumerable<int[]> GetVariableModificationPatterns(List<KeyValuePair<int, List<ModificationWithMass>>> possibleVariableModifications, int unmodifiedResiduesDesired, int[] variableModificationPattern, int index)
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

        private static Dictionary<int, ModificationWithMass> GetNewVariableModificationPattern(int[] variableModificationArray, IEnumerable<KeyValuePair<int, List<ModificationWithMass>>> possibleVariableModifications)
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

        #endregion Private Methods
    }
}