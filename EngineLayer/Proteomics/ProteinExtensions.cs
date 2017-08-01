using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public static class ProteinExtensions
    {

        #region Public Methods

        public static IEnumerable<PeptideWithPossibleModifications> Digest(this Protein protein, Protease protease, int maximumMissedCleavages, int? minPeptidesLength, int? maxPeptidesLength, InitiatorMethionineBehavior initiatorMethionineBehavior, IEnumerable<ModificationWithMass> allKnownFixedModifications, bool addCompIons)
        {
            if (protease.CleavageSpecificity != CleavageSpecificity.None)
            {
                // these are the 1-based residue indices the protease cleaves AFTER
                List<int> oneBasedIndicesToCleaveAfter = protease.GetDigestionSiteIndices(protein.BaseSequence);
                // Cleave after 0, or before the first one
                oneBasedIndicesToCleaveAfter.Insert(0, 0);
                // Cleave after Length index
                oneBasedIndicesToCleaveAfter.Add(protein.Length);

                if (protease.CleavageSpecificity == CleavageSpecificity.Full)
                {
                    for (int missed_cleavages = 0; missed_cleavages <= maximumMissedCleavages; missed_cleavages++)
                    {
                        for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - missed_cleavages - 1; i++)
                        {
                            // Retain!
                            if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave || i != 0 || protein[0] != 'M')
                            {
                                if ((!minPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1] - oneBasedIndicesToCleaveAfter[i] >= minPeptidesLength) &&
                                    (!maxPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1] - oneBasedIndicesToCleaveAfter[i] <= maxPeptidesLength))
                                {
                                    yield return new PeptideWithPossibleModifications(oneBasedIndicesToCleaveAfter[i] + 1, oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1], protein, missed_cleavages, "full", allKnownFixedModifications, addCompIons);
                                }
                            }
                            // Cleave!
                            if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && i == 0 && protein[0] == 'M')
                            {
                                if ((!minPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1] - 1 >= minPeptidesLength) &&
                                    (!maxPeptidesLength.HasValue || oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1] - 1 <= maxPeptidesLength))
                                {
                                    yield return new PeptideWithPossibleModifications(2, oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1], protein, missed_cleavages, "full:M cleaved", allKnownFixedModifications, addCompIons);
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
                                        yield return new PeptideWithPossibleModifications(proteolysisProduct.OneBasedBeginPosition.Value, oneBasedIndicesToCleaveAfter[i + missed_cleavages], protein, missed_cleavages, proteolysisProduct.Type + " start", allKnownFixedModifications, addCompIons);
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
                                        yield return new PeptideWithPossibleModifications(oneBasedIndicesToCleaveAfter[i - missed_cleavages - 1] + 1, proteolysisProduct.OneBasedEndPosition.Value, protein, missed_cleavages, proteolysisProduct.Type + " end", allKnownFixedModifications, addCompIons);
                                    }
                                }
                            }
                    }
                }
                else // protease.CleavageSpecificity == CleavageSpecificity.Semi
                {
                    throw new NotImplementedException();
                }
            }
            else if (protease.Name.Equals("singleC"))
            {
                //cleave in one spot
                for (int index = 1; index <= protein.Length; index++) //position BEFORE the amino acid
                {
                    if ((!minPeptidesLength.HasValue || (index) >= minPeptidesLength))//&&
                    {
                        yield return new PeptideWithPossibleModifications(Math.Max(1, index - (maxPeptidesLength ?? 50)), index, protein, 0, "MHC", allKnownFixedModifications, addCompIons);
                    }
                }
            }
            else if (protease.Name.Equals("singleN"))
            {
                //cleave in one spot
                for (int index = 1; index <= protein.Length; index++) //position BEFORE the amino acid
                {
                    if ((!minPeptidesLength.HasValue || (protein.Length - index + 1) >= minPeptidesLength)) //&&
                    {
                        yield return new PeptideWithPossibleModifications(index, Math.Min(protein.Length, index + (maxPeptidesLength ?? 50)), protein, 0, "MHC", allKnownFixedModifications, addCompIons);
                    }
                }
            }
            else  // protease.CleavageSpecificity == CleavageSpecificity.None
            {
                if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave || protein[0] != 'M')
                {
                    if ((!minPeptidesLength.HasValue || protein.Length >= minPeptidesLength) &&
                                    (!maxPeptidesLength.HasValue || protein.Length <= maxPeptidesLength))
                    {
                        yield return new PeptideWithPossibleModifications(1, protein.Length, protein, 0, "full", allKnownFixedModifications, addCompIons);
                    }
                }
                if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && protein[0] == 'M')
                {
                    if ((!minPeptidesLength.HasValue || protein.Length - 1 >= minPeptidesLength) &&
                                    (!maxPeptidesLength.HasValue || protein.Length - 1 <= maxPeptidesLength))
                    {
                        yield return new PeptideWithPossibleModifications(2, protein.Length, protein, 0, "full:M cleaved", allKnownFixedModifications, addCompIons);
                    }
                }

                // Also digest using the proteolysis product start/end indices
                foreach (var proteolysisProduct in protein.ProteolysisProducts)
                    if (proteolysisProduct.OneBasedEndPosition.HasValue && proteolysisProduct.OneBasedBeginPosition.HasValue)
                    {
                        if ((!minPeptidesLength.HasValue || (proteolysisProduct.OneBasedEndPosition.Value - proteolysisProduct.OneBasedBeginPosition.Value + 1) >= minPeptidesLength) &&
                                       (!maxPeptidesLength.HasValue || (proteolysisProduct.OneBasedEndPosition.Value - proteolysisProduct.OneBasedBeginPosition.Value + 1) <= maxPeptidesLength))
                        {
                            yield return new PeptideWithPossibleModifications(proteolysisProduct.OneBasedBeginPosition.Value, proteolysisProduct.OneBasedEndPosition.Value, protein, 0, proteolysisProduct.Type, allKnownFixedModifications, addCompIons);
                        }
                    }
            }
        }

        public static string TabSeparatedString(this Protein protein)
        {
            return protein.Accession + "\t" + protein.FullName + "\t" + string.Join(",", protein.GeneNames.Select(c => c.Item1 + ":" + c.Item2));
        }

        #endregion Public Methods

    }
}