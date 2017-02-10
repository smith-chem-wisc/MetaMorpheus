using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public static class ProteinExtensions
    {

        public static IEnumerable<PeptideWithPossibleModifications> Digest(this Protein protein, Protease protease, int maximumMissedCleavages, InitiatorMethionineBehavior initiatorMethionineBehavior, IEnumerable<ModificationWithMass> allKnownFixedModifications)
        {
            var k = GetListWithMass(protein.OneBasedPossibleLocalizedModifications).ToList();
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
                                yield return new PeptideWithPossibleModifications(oneBasedIndicesToCleaveAfter[i] + 1, oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1], protein, missed_cleavages, "full", allKnownFixedModifications, k);
                            }
                            // Cleave!
                            if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && i == 0 && protein[0] == 'M')
                            {
                                yield return new PeptideWithPossibleModifications(2, oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1], protein, missed_cleavages, "full:M cleaved", allKnownFixedModifications, k);
                            }
                        }

                        // Also digest using the chain peptide start/end indices
                        for (int chainPeptideIndex = 0; chainPeptideIndex < protein.OneBasedBeginPositions.Length; chainPeptideIndex++)
                        {
                            if (protein.OneBasedBeginPositions[chainPeptideIndex] != 1 || protein.OneBasedEndPositions[chainPeptideIndex] != protein.Length)
                            {
                                int i = 0;
                                while (oneBasedIndicesToCleaveAfter[i] < protein.OneBasedBeginPositions[chainPeptideIndex])
                                    i++;
                                // Start peptide
                                if (i + missed_cleavages < oneBasedIndicesToCleaveAfter.Count && oneBasedIndicesToCleaveAfter[i + missed_cleavages] <= protein.OneBasedEndPositions[chainPeptideIndex] && protein.OneBasedBeginPositions[chainPeptideIndex].HasValue)
                                    yield return new PeptideWithPossibleModifications(protein.OneBasedBeginPositions[chainPeptideIndex].Value, oneBasedIndicesToCleaveAfter[i + missed_cleavages], protein, missed_cleavages, protein.BigPeptideTypes[chainPeptideIndex] + " start", allKnownFixedModifications, k);

                                while (oneBasedIndicesToCleaveAfter[i] < protein.OneBasedEndPositions[chainPeptideIndex])
                                    i++;
                                // End
                                if (i - missed_cleavages - 1 >= 0 && oneBasedIndicesToCleaveAfter[i - missed_cleavages - 1] + 1 >= protein.OneBasedBeginPositions[chainPeptideIndex] && protein.OneBasedEndPositions[chainPeptideIndex].HasValue)
                                    yield return new PeptideWithPossibleModifications(oneBasedIndicesToCleaveAfter[i - missed_cleavages - 1] + 1, protein.OneBasedEndPositions[chainPeptideIndex].Value, protein, missed_cleavages, protein.BigPeptideTypes[chainPeptideIndex] + " end", allKnownFixedModifications, k);
                            }
                        }
                    }
                }
                else
                {
                    throw new NotImplementedException();
                }
            }
            else  // protease.CleavageSpecificity == CleavageSpecificity.None
            {
                if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave || protein[0] != 'M')
                {
                    yield return new PeptideWithPossibleModifications(1, protein.Length, protein, 0, "full", allKnownFixedModifications, k);
                }
                if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && protein[0] == 'M')
                {
                    yield return new PeptideWithPossibleModifications(2, protein.Length, protein, 0, "full:M cleaved", allKnownFixedModifications, k);
                }

                // Also digest using the chain peptide start/end indices
                for (int chainPeptideIndex = 0; chainPeptideIndex < protein.OneBasedBeginPositions.Length; chainPeptideIndex++)
                {
                    if (protein.OneBasedEndPositions[chainPeptideIndex].HasValue && protein.OneBasedBeginPositions[chainPeptideIndex].HasValue)
                        yield return new PeptideWithPossibleModifications(protein.OneBasedBeginPositions[chainPeptideIndex].Value, protein.OneBasedEndPositions[chainPeptideIndex].Value, protein, 0, protein.BigPeptideTypes[chainPeptideIndex], allKnownFixedModifications, k);
                }
            }
        }

        private static IEnumerable<Tuple<int, List<ModificationWithMass>>> GetListWithMass(IDictionary<int, List<Modification>> oneBasedPossibleLocalizedModifications)
        {
            foreach (var hm in oneBasedPossibleLocalizedModifications)
            {
                List<ModificationWithMass> list = new List<ModificationWithMass>();
                foreach (var ok in hm.Value)
                {
                    if (ok is ModificationWithMass)
                        list.Add(ok as ModificationWithMass);
                }
                yield return new Tuple<int, List<ModificationWithMass>>(hm.Key, list);
            }
        }
    }
}