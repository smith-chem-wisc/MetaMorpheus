using Proteomics;
using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public static class ProteinExtensions
    {

        #region Public Methods

        public static IEnumerable<PeptideWithPossibleModifications> Digest(this Protein protein, Protease protease, int maximumMissedCleavages, int? minPeptidesLenght, int? maxPeptidesLength, InitiatorMethionineBehavior initiatorMethionineBehavior, IEnumerable<ModificationWithMass> allKnownFixedModifications)
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
                                PeptideWithPossibleModifications pwpm = new PeptideWithPossibleModifications(oneBasedIndicesToCleaveAfter[i] + 1, oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1], protein, missed_cleavages, "full", allKnownFixedModifications);
                                if (!minPeptidesLenght.HasValue && !maxPeptidesLength.HasValue) { yield return pwpm; }
                                else if (pwpm.Length >= minPeptidesLenght && !maxPeptidesLength.HasValue) { yield return pwpm; }
                                else if (pwpm.Length >= minPeptidesLenght && pwpm.Length <= maxPeptidesLength) { yield return pwpm; }
                                else if (!minPeptidesLenght.HasValue && pwpm.Length <= maxPeptidesLength) { yield return pwpm; }
                            }
                            // Cleave!
                            if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && i == 0 && protein[0] == 'M')
                            {
                                PeptideWithPossibleModifications pwpm = new PeptideWithPossibleModifications(2, oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1], protein, missed_cleavages, "full:M cleaved", allKnownFixedModifications);
                                if (!minPeptidesLenght.HasValue && !maxPeptidesLength.HasValue) { yield return pwpm; }
                                else if (pwpm.Length >= minPeptidesLenght && !maxPeptidesLength.HasValue) { yield return pwpm; }
                                else if (pwpm.Length >= minPeptidesLenght && pwpm.Length <= maxPeptidesLength) { yield return pwpm; }
                                else if (!minPeptidesLenght.HasValue && pwpm.Length <= maxPeptidesLength) { yield return pwpm; }
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
                                    var pwpm = new PeptideWithPossibleModifications(proteolysisProduct.OneBasedBeginPosition.Value, oneBasedIndicesToCleaveAfter[i + missed_cleavages], protein, missed_cleavages, proteolysisProduct.Type + " start", allKnownFixedModifications);
                                    if (!minPeptidesLenght.HasValue && !maxPeptidesLength.HasValue) { yield return pwpm; }
                                    else if (pwpm.Length >= minPeptidesLenght && !maxPeptidesLength.HasValue) { yield return pwpm; }
                                    else if (pwpm.Length >= minPeptidesLenght && pwpm.Length <= maxPeptidesLength) { yield return pwpm; }
                                    else if (!minPeptidesLenght.HasValue && pwpm.Length <= maxPeptidesLength) { yield return pwpm; }
                                }
                                while (oneBasedIndicesToCleaveAfter[i] < proteolysisProduct.OneBasedEndPosition)
                                    i++;
                                // End
                                if (i - missed_cleavages - 1 >= 0 && oneBasedIndicesToCleaveAfter[i - missed_cleavages - 1] + 1 >= proteolysisProduct.OneBasedBeginPosition && proteolysisProduct.OneBasedEndPosition.HasValue)
                                {
                                    var pwpm = new PeptideWithPossibleModifications(oneBasedIndicesToCleaveAfter[i - missed_cleavages - 1] + 1, proteolysisProduct.OneBasedEndPosition.Value, protein, missed_cleavages, proteolysisProduct.Type + " end", allKnownFixedModifications);
                                    if (!minPeptidesLenght.HasValue && !maxPeptidesLength.HasValue) { yield return pwpm; }
                                    else if (pwpm.Length >= minPeptidesLenght && !maxPeptidesLength.HasValue) { yield return pwpm; }
                                    else if (pwpm.Length >= minPeptidesLenght && pwpm.Length <= maxPeptidesLength) { yield return pwpm; }
                                    else if (!minPeptidesLenght.HasValue && pwpm.Length <= maxPeptidesLength) { yield return pwpm; }
                                }
                            }
                    }
                }
                else // protease.CleavageSpecificity == CleavageSpecificity.Semi
                {
                    throw new NotImplementedException();
                }
            }
            else  // protease.CleavageSpecificity == CleavageSpecificity.None
            {
                if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave || protein[0] != 'M')
                {
                    var pwpm = new PeptideWithPossibleModifications(1, protein.Length, protein, 0, "full", allKnownFixedModifications);
                    if (!minPeptidesLenght.HasValue && !maxPeptidesLength.HasValue) { yield return pwpm; }
                    else if (pwpm.Length >= minPeptidesLenght && !maxPeptidesLength.HasValue) { yield return pwpm; }
                    else if (pwpm.Length >= minPeptidesLenght && pwpm.Length <= maxPeptidesLength) { yield return pwpm; }
                    else if (!minPeptidesLenght.HasValue && pwpm.Length <= maxPeptidesLength) { yield return pwpm; }
                }
                if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && protein[0] == 'M')
                {
                    var pwpm = new PeptideWithPossibleModifications(2, protein.Length, protein, 0, "full:M cleaved", allKnownFixedModifications);
                    if (!minPeptidesLenght.HasValue && !maxPeptidesLength.HasValue) { yield return pwpm; }
                    else if (pwpm.Length >= minPeptidesLenght && !maxPeptidesLength.HasValue) { yield return pwpm; }
                    else if (pwpm.Length >= minPeptidesLenght && pwpm.Length <= maxPeptidesLength) { yield return pwpm; }
                    else if (!minPeptidesLenght.HasValue && pwpm.Length <= maxPeptidesLength) { yield return pwpm; }
                }

                // Also digest using the proteolysis product start/end indices
                foreach (var proteolysisProduct in protein.ProteolysisProducts)
                    if (proteolysisProduct.OneBasedEndPosition.HasValue && proteolysisProduct.OneBasedBeginPosition.HasValue)
                    {
                        var pwpm = new PeptideWithPossibleModifications(proteolysisProduct.OneBasedBeginPosition.Value, proteolysisProduct.OneBasedEndPosition.Value, protein, 0, proteolysisProduct.Type, allKnownFixedModifications);
                        if (!minPeptidesLenght.HasValue && !maxPeptidesLength.HasValue) { yield return pwpm; }
                        else if (pwpm.Length >= minPeptidesLenght && !maxPeptidesLength.HasValue) { yield return pwpm; }
                        else if (pwpm.Length >= minPeptidesLenght && pwpm.Length <= maxPeptidesLength) { yield return pwpm; }
                        else if (!minPeptidesLenght.HasValue && pwpm.Length <= maxPeptidesLength) { yield return pwpm; }
                    }
            }
        }
        #endregion Public Methods

    }
}