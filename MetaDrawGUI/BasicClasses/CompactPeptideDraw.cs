using System;
using System.Linq;
using EngineLayer;
using Proteomics;
using System.Collections.Generic;

namespace MetaDrawGUI
{

    [Serializable]
    public class CompactPeptideDraw : CompactPeptideBase
    {
        private const int digitsForRoundingMasses = 7;

        public CompactPeptideDraw(PepWithSetModForCompactPep peptideWithSetModifications, TerminusType terminusType)
        {
            NTerminalMasses = null;
            CTerminalMasses = null;
            if (terminusType == TerminusType.None || terminusType == TerminusType.N)
            {
                NTerminalMasses = ComputeFollowingFragmentMassesByPepCom(peptideWithSetModifications, 0, 0, 1).ToArray();
            }
            if (terminusType == TerminusType.None || terminusType == TerminusType.C)
            {
                CTerminalMasses = ComputeFollowingFragmentMassesByPepCom(peptideWithSetModifications, 0, peptideWithSetModifications.Length + 1, -1).ToArray();
            }
            MonoisotopicMassIncludingFixedMods = peptideWithSetModifications.MonoisotopicMass;
        }

        protected static IEnumerable<double> ComputeFollowingFragmentMassesByPepCom(PepWithSetModForCompactPep yyy, double prevMass, int oneBasedIndexToLookAt, int direction)
        {
            ModificationWithMass currentModification = null;
            do
            {
                if (oneBasedIndexToLookAt != 0 && oneBasedIndexToLookAt != yyy.Length + 1)
                {
                    prevMass += Residue.ResidueMonoisotopicMass[yyy.BaseSequence[oneBasedIndexToLookAt - 1]];
                }

                // If modification exists
                if (yyy.allModsOneIsNterminus.TryGetValue(oneBasedIndexToLookAt + 1, out currentModification))
                {
                    if (currentModification.neutralLosses.Count == 1 && oneBasedIndexToLookAt != 0 && oneBasedIndexToLookAt != yyy.Length + 1)
                    {
                        prevMass += currentModification.monoisotopicMass - currentModification.neutralLosses.First();
                        yield return Math.Round(prevMass, digitsForRoundingMasses);
                    }
                    else
                    {
                        foreach (double nl in currentModification.neutralLosses)
                        {
                            var theMass = prevMass + currentModification.monoisotopicMass - nl;
                            if (oneBasedIndexToLookAt != 0 && oneBasedIndexToLookAt != yyy.Length + 1)
                            {
                                yield return Math.Round(theMass, digitsForRoundingMasses);
                            }
                            if ((direction == 1 && oneBasedIndexToLookAt + direction < yyy.Length) ||
                                (direction == -1 && oneBasedIndexToLookAt + direction > 1))
                            {
                                foreach (var nextMass in ComputeFollowingFragmentMassesByPepCom(yyy, theMass, oneBasedIndexToLookAt + direction, direction))
                                {
                                    yield return Math.Round(nextMass, digitsForRoundingMasses);
                                }
                            }
                        }
                        break;
                    }
                }
                else if (oneBasedIndexToLookAt != 0 && oneBasedIndexToLookAt != yyy.Length + 1) // No modification exists
                {
                    yield return Math.Round(prevMass, digitsForRoundingMasses);
                }
                oneBasedIndexToLookAt += direction;
            } while ((oneBasedIndexToLookAt > 1 && direction == -1) || (oneBasedIndexToLookAt < yyy.Length && direction == 1));
        }
    }
}