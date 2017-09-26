using Chemistry;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    [Serializable]
    public abstract class CompactPeptideBase : IEquatable<CompactPeptideBase>
    {
        #region Protected Fields

        protected static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;

        protected static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        protected static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;

        protected static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        #endregion Protected Fields

        #region Private Fields

        private const int digitsForRoundingMasses = 7;

        private const double massTolForPeptideEquality = 1e-7;

        #endregion Private Fields

        #region Public Properties

        public double[] CTerminalMasses { get; protected set; }
        public double[] NTerminalMasses { get; protected set; }
        public double MonoisotopicMassIncludingFixedMods { get; protected set; }

        #endregion Public Properties

        #region Public Methods

        public double[] ProductMassesMightHaveDuplicatesAndNaNs(List<ProductType> productTypes)
        {
            int massLen = 0;
            bool containsAdot = productTypes.Contains(ProductType.Adot);
            bool containsB = productTypes.Contains(ProductType.B);
            bool containsBnoB1 = productTypes.Contains(ProductType.BnoB1ions);
            bool containsC = productTypes.Contains(ProductType.C);
            bool containsX = productTypes.Contains(ProductType.X);
            bool containsY = productTypes.Contains(ProductType.Y);
            bool containsZdot = productTypes.Contains(ProductType.Zdot);

            if (containsAdot)
                throw new NotImplementedException();
            if (containsBnoB1)
                massLen += NTerminalMasses.Length - 1;
            else if (containsB)
                massLen += NTerminalMasses.Length;
            if (containsC)
                massLen += NTerminalMasses.Length;
            if (containsX)
                throw new NotImplementedException();
            if (containsY)
                massLen += CTerminalMasses.Length;
            if (containsZdot)
                massLen += CTerminalMasses.Length;

            double[] massesToReturn = new double[massLen];

            int i = 0;
            if (NTerminalMasses != null)
                for (int j = 0; j < NTerminalMasses.Length; j++)
                {
                    var hm = NTerminalMasses[j];
                    if (containsBnoB1 && j > 0)
                    {
                        massesToReturn[i] = hm;
                        i++;
                    }
                    else if (containsB)
                    {
                        massesToReturn[i] = hm;
                        i++;
                    }
                    if (containsC)
                    {
                        massesToReturn[i] = hm + nitrogenAtomMonoisotopicMass + 3 * hydrogenAtomMonoisotopicMass;
                        i++;
                    }
                }
            if (CTerminalMasses != null)
                foreach (double hm in CTerminalMasses)
                {
                    if (containsY)
                    {
                        massesToReturn[i] = hm + waterMonoisotopicMass;
                        i++;
                    }
                    if (containsZdot)
                    {
                        massesToReturn[i] = hm + oxygenAtomMonoisotopicMass - nitrogenAtomMonoisotopicMass;
                        i++;
                    }
                }
            return massesToReturn;
        }

        /// <summary>
        /// Sometimes says not equal when in reality should be equal, due to rounding errors. Small but annoying bug. Careful when fixing! Make sure Indexing runs at a reasonable speed.
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public override bool Equals(object obj)
        {
            var cp = obj as CompactPeptideBase;
            return cp != null && Equals(cp);
        }

        public override int GetHashCode()
        {
            unchecked
            {
                {
                    var result = 0;
                    if (CTerminalMasses == null)
                    {
                        foreach (double b in NTerminalMasses)
                            result = (result * 31) ^ b.GetHashCode();
                    }
                    else
                    {
                        foreach (double b in CTerminalMasses)
                            result = (result * 31) ^ b.GetHashCode();
                    }
                    return result;
                }
            }
        }

        public bool Equals(CompactPeptideBase cp)
        {
            if (CTerminalMasses != null && cp.CTerminalMasses != null)
            {
                if (NTerminalMasses != null && cp.NTerminalMasses != null) //neither series is nulll
                {
                    return (
                        ((double.IsNaN(MonoisotopicMassIncludingFixedMods) && double.IsNaN(cp.MonoisotopicMassIncludingFixedMods)) || (Math.Abs(MonoisotopicMassIncludingFixedMods - cp.MonoisotopicMassIncludingFixedMods) < massTolForPeptideEquality))
                        && CTerminalMasses.SequenceEqual(cp.CTerminalMasses)
                        && NTerminalMasses.SequenceEqual(cp.NTerminalMasses)
                        );
                }
                else //No N-terminal ions
                {
                    return (
                        ((double.IsNaN(MonoisotopicMassIncludingFixedMods) && double.IsNaN(cp.MonoisotopicMassIncludingFixedMods)) || (Math.Abs(MonoisotopicMassIncludingFixedMods - cp.MonoisotopicMassIncludingFixedMods) < massTolForPeptideEquality))
                        && CTerminalMasses.SequenceEqual(cp.CTerminalMasses)
                        );
                }
            }
            else if (NTerminalMasses != null && cp.NTerminalMasses != null) //No C-terminal ions
            {
                return (
                    ((double.IsNaN(MonoisotopicMassIncludingFixedMods) && double.IsNaN(cp.MonoisotopicMassIncludingFixedMods)) || (Math.Abs(MonoisotopicMassIncludingFixedMods - cp.MonoisotopicMassIncludingFixedMods) < massTolForPeptideEquality))
                    && NTerminalMasses.SequenceEqual(cp.NTerminalMasses)
                    );
            }
            else //Cannot compare
            {
                return false;
            }
        }

        #endregion Public Methods

        #region Protected Methods

        protected static IEnumerable<double> ComputeFollowingFragmentMasses(PeptideWithSetModifications yyy, double prevMass, int oneBasedIndexToLookAt, int direction)
        {
            ModificationWithMass currentModification = null;
            do
            {
                if (oneBasedIndexToLookAt != 0 && oneBasedIndexToLookAt != yyy.Length + 1)
                    prevMass += Residue.ResidueMonoisotopicMass[yyy[oneBasedIndexToLookAt - 1]];

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
                                yield return Math.Round(theMass, digitsForRoundingMasses);
                            if ((direction == 1 && oneBasedIndexToLookAt + direction < yyy.Length) ||
                                (direction == -1 && oneBasedIndexToLookAt + direction > 1))
                                foreach (var nextMass in ComputeFollowingFragmentMasses(yyy, theMass, oneBasedIndexToLookAt + direction, direction))
                                    yield return Math.Round(nextMass, digitsForRoundingMasses);
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

        #endregion Protected Methods
    }
}