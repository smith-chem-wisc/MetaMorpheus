using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public static class ProteinExtensions
    {
        #region Public Methods

        public static IEnumerable<PeptideWithPossibleModifications> Digest(this Protein protein, DigestionParams digestionParams, IEnumerable<ModificationWithMass> allKnownFixedModifications)
        {
            var protease = digestionParams.Protease;
            var maximumMissedCleavages = digestionParams.MaxMissedCleavages;
            var initiatorMethionineBehavior = digestionParams.InitiatorMethionineBehavior;
            var minPeptidesLength = digestionParams.MinPeptideLength;
            var maxPeptidesLength = digestionParams.MaxPeptideLength;
            List<int> oneBasedIndicesToCleaveAfter = new List<int>();

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
                                    yield return new PeptideWithPossibleModifications(oneBasedIndicesToCleaveAfter[i] + 1, oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1], protein, missed_cleavages, "full", allKnownFixedModifications);
                                }
                            }
                            // Cleave!
                            if (i == 0 && initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && protein[0] == 'M')
                            {
                                if ((!minPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1] - 1 >= minPeptidesLength) &&
                                    (!maxPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1] - 1 <= maxPeptidesLength))
                                {
                                    yield return new PeptideWithPossibleModifications(2, oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1], protein, missed_cleavages, "full:M cleaved", allKnownFixedModifications);
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
                                        yield return new PeptideWithPossibleModifications(proteolysisProduct.OneBasedBeginPosition.Value, oneBasedIndicesToCleaveAfter[i + missed_cleavages], protein, missed_cleavages, proteolysisProduct.Type + " start", allKnownFixedModifications);
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
                                        yield return new PeptideWithPossibleModifications(oneBasedIndicesToCleaveAfter[i - missed_cleavages - 1] + 1, proteolysisProduct.OneBasedEndPosition.Value, protein, missed_cleavages, proteolysisProduct.Type + " end", allKnownFixedModifications);
                                    }
                                }
                            }
                    }
                    break;

                case CleavageSpecificity.FullMaxN:
                case CleavageSpecificity.FullMaxC:
                    // these are the 1-based residue indices the protease cleaves AFTER
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
                            yield return new PeptideWithPossibleModifications(oneBasedIndicesToCleaveAfter[i] + 1, oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1], protein, oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1] - oneBasedIndicesToCleaveAfter[i], "semi", allKnownFixedModifications);

                        // Cleave!
                        if ((i == 0 && initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && protein[0] == 'M')
                            && (!minPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1] - 1 >= minPeptidesLength)
                            && (!maxPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1] - 1 <= maxPeptidesLength))
                            yield return new PeptideWithPossibleModifications(2, oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1], protein, oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1] - 1, "semi:M cleaved", allKnownFixedModifications);
                    }

                    int lastIndex = oneBasedIndicesToCleaveAfter.Count - 1;
                    int maxIndex = maximumMissedCleavages < lastIndex ? maximumMissedCleavages : lastIndex;
                    if (protease.CleavageSpecificity == CleavageSpecificity.FullMaxN) //tricky, it's N because we want the extra peptide at the C terminus |_
                    {
                        for (int i = 1; i <= maxIndex; i++)
                        {
                            if ((!minPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[lastIndex] - oneBasedIndicesToCleaveAfter[lastIndex - i] >= minPeptidesLength) &&
                                (!maxPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[lastIndex] - oneBasedIndicesToCleaveAfter[lastIndex - i] <= maxPeptidesLength))
                            {
                                yield return new PeptideWithPossibleModifications(oneBasedIndicesToCleaveAfter[lastIndex - i] + 1, oneBasedIndicesToCleaveAfter[lastIndex], protein, oneBasedIndicesToCleaveAfter[lastIndex] - oneBasedIndicesToCleaveAfter[lastIndex - i], "semi", allKnownFixedModifications);
                            }
                        }
                    }
                    else //FullMaxC, never cleave M
                    {
                        for (int i = 1; i <= maxIndex; i++)
                        {
                            if ((!minPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i] - oneBasedIndicesToCleaveAfter[0] >= minPeptidesLength) &&
                                (!maxPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i] - oneBasedIndicesToCleaveAfter[0] <= maxPeptidesLength))
                            {
                                yield return new PeptideWithPossibleModifications(oneBasedIndicesToCleaveAfter[0] + 1, oneBasedIndicesToCleaveAfter[i], protein, oneBasedIndicesToCleaveAfter[i] - oneBasedIndicesToCleaveAfter[0], "semi", allKnownFixedModifications);
                            }
                        }
                    }
                    // Also digest using the proteolysis product start/end indices
                    foreach (var proteolysisProduct in protein.ProteolysisProducts)
                        if (proteolysisProduct.OneBasedBeginPosition != 1 || proteolysisProduct.OneBasedEndPosition != protein.Length)
                        {
                            yield return new PeptideWithPossibleModifications(proteolysisProduct.OneBasedBeginPosition.Value, proteolysisProduct.OneBasedEndPosition.Value, protein, 0, proteolysisProduct.Type + " start", allKnownFixedModifications);
                        }
                    break;

                case CleavageSpecificity.SingleN:
                    //cleave in one spot
                    for (int index = 1; index <= protein.Length; index++) //position BEFORE the amino acid
                    {
                        if ((!minPeptidesLength.HasValue || (protein.Length - index + 1) >= minPeptidesLength)) //&&
                        {
                            yield return new PeptideWithPossibleModifications(index, Math.Min(protein.Length, index + (maxPeptidesLength ?? 50)), protein, 0, "SingleN", allKnownFixedModifications);
                        }
                    }
                    break;

                case CleavageSpecificity.SingleC:
                    //cleave in one spot
                    for (int index = 1; index <= protein.Length; index++) //position BEFORE the amino acid
                    {
                        if ((!minPeptidesLength.HasValue || (index) >= minPeptidesLength))//&&
                        {
                            yield return new PeptideWithPossibleModifications(Math.Max(1, index - (maxPeptidesLength ?? 50)), index, protein, 0, "SingleC", allKnownFixedModifications);
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
                                yield return new PeptideWithPossibleModifications(nTermIndexofProtein + 1, cTermIndexofProtein, protein, cTermIndexofProtein - nTermIndexofProtein, "semi", allKnownFixedModifications);
                            }
                            //fixedC
                            for (int j = nTermIndexofProtein + 1; j < cTermIndexofProtein; j++)
                            {
                                if ((!minPeptidesLength.HasValue || cTermIndexofProtein - j >= minPeptidesLength) &&
                                    (!maxPeptidesLength.HasValue || cTermIndexofProtein - j <= maxPeptidesLength))
                                {
                                    yield return new PeptideWithPossibleModifications(j + 1, cTermIndexofProtein, protein, cTermIndexofProtein - j, "semiC", allKnownFixedModifications);
                                }
                            }
                            //fixedN
                            for (int j = cTermIndexofProtein - 1; j > nTermIndexofProtein; j--)
                            {
                                if ((!minPeptidesLength.HasValue || j - nTermIndexofProtein >= minPeptidesLength) &&
                                    (!maxPeptidesLength.HasValue || j - nTermIndexofProtein <= maxPeptidesLength))
                                {
                                    yield return new PeptideWithPossibleModifications(nTermIndexofProtein + 1, j, protein, j - nTermIndexofProtein, "semiN", allKnownFixedModifications);
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
                                yield return new PeptideWithPossibleModifications(nTermIndexOfProtein + 1, cTermIndexOfProtein, protein, cTermIndexOfProtein - nTermIndexOfProtein, "semi:M cleaved", allKnownFixedModifications);
                            }
                            //fixedC
                            for (int j = nTermIndexOfProtein + 1; j < cTermIndexOfProtein; j++)
                            {
                                if ((!minPeptidesLength.HasValue || cTermIndexOfProtein - j >= minPeptidesLength) &&
                                    (!maxPeptidesLength.HasValue || cTermIndexOfProtein - j <= maxPeptidesLength))
                                {
                                    yield return new PeptideWithPossibleModifications(j + 1, cTermIndexOfProtein, protein, cTermIndexOfProtein - j, "semiC", allKnownFixedModifications);
                                }
                            }
                            //fixedN
                            for (int j = cTermIndexOfProtein - 1; j > nTermIndexOfProtein; j--)
                            {
                                if ((!minPeptidesLength.HasValue || j - nTermIndexOfProtein >= minPeptidesLength) &&
                                    (!maxPeptidesLength.HasValue || j - nTermIndexOfProtein <= maxPeptidesLength))
                                {
                                    yield return new PeptideWithPossibleModifications(nTermIndexOfProtein + 1, j, protein, j - nTermIndexOfProtein, "semiN", allKnownFixedModifications);
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
                                yield return new PeptideWithPossibleModifications(nTermIndexOfProtein + 1, j, protein, j - nTermIndexOfProtein, "semiN", allKnownFixedModifications);
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
                                yield return new PeptideWithPossibleModifications(j + 1, cTermIndexOfProtein, protein, cTermIndexOfProtein - j, "semiC", allKnownFixedModifications);
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
                                    yield return new PeptideWithPossibleModifications(proteolysisProduct.OneBasedBeginPosition.Value, j, protein, j - proteolysisProduct.OneBasedBeginPosition.Value, proteolysisProduct.Type + " start", allKnownFixedModifications);
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
                                    yield return new PeptideWithPossibleModifications(j, proteolysisProduct.OneBasedEndPosition.Value, protein, proteolysisProduct.OneBasedEndPosition.Value - j, proteolysisProduct.Type + " end", allKnownFixedModifications);
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
                            yield return new PeptideWithPossibleModifications(1, protein.Length, protein, 0, "full", allKnownFixedModifications);
                        }
                    }
                    if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && protein[0] == 'M')
                    {
                        if ((!minPeptidesLength.HasValue || protein.Length - 1 >= minPeptidesLength) &&
                                        (!maxPeptidesLength.HasValue || protein.Length - 1 <= maxPeptidesLength))
                        {
                            yield return new PeptideWithPossibleModifications(2, protein.Length, protein, 0, "full:M cleaved", allKnownFixedModifications);
                        }
                    }

                    // Also digest using the proteolysis product start/end indices
                    foreach (var proteolysisProduct in protein.ProteolysisProducts)
                        if (proteolysisProduct.OneBasedEndPosition.HasValue && proteolysisProduct.OneBasedBeginPosition.HasValue)
                        {
                            if ((!minPeptidesLength.HasValue || (proteolysisProduct.OneBasedEndPosition.Value - proteolysisProduct.OneBasedBeginPosition.Value + 1) >= minPeptidesLength) &&
                                           (!maxPeptidesLength.HasValue || (proteolysisProduct.OneBasedEndPosition.Value - proteolysisProduct.OneBasedBeginPosition.Value + 1) <= maxPeptidesLength))
                            {
                                yield return new PeptideWithPossibleModifications(proteolysisProduct.OneBasedBeginPosition.Value, proteolysisProduct.OneBasedEndPosition.Value, protein, 0, proteolysisProduct.Type, allKnownFixedModifications);
                            }
                        }
                    break;

                default:
                    throw new NotImplementedException();
            }
        }

        public static string TabSeparatedString(this Protein protein)
        {
            return protein.Accession + "\t" + protein.FullName + "\t" + string.Join(",", protein.GeneNames.Select(c => c.Item1 + ":" + c.Item2));
        }

        #endregion Public Methods
    }
}