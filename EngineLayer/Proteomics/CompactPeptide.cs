using Chemistry;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    [Serializable]
    public class CompactPeptide
    {

        #region Private Fields

        private static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        private static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;
        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        #endregion Private Fields

        #region Public Constructors

        public CompactPeptide(PeptideWithSetModifications peptideWithSetModifications)
        {
            ModificationWithMass pep_n_term_variable_mod;
            double theMass = 0;
            if (peptideWithSetModifications.allModsOneIsNterminus.TryGetValue(1, out pep_n_term_variable_mod))
                foreach (double nl in pep_n_term_variable_mod.neutralLosses)
                    theMass = pep_n_term_variable_mod.monoisotopicMass - nl;
            else
                theMass = 0;
            NTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, theMass, 1, 1).ToArray();

            ModificationWithMass pep_c_term_variable_mod;
            theMass = 0;
            if (peptideWithSetModifications.allModsOneIsNterminus.TryGetValue(peptideWithSetModifications.Length + 2, out pep_c_term_variable_mod))
                foreach (double nl in pep_c_term_variable_mod.neutralLosses)
                    theMass = pep_c_term_variable_mod.monoisotopicMass - nl;
            else
                theMass = 0;
            CTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, theMass, peptideWithSetModifications.Length, -1).ToArray();

            MonoisotopicMassIncludingFixedMods = peptideWithSetModifications.MonoisotopicMass;
        }

        #endregion Public Constructors

        #region Public Properties

        public double[] CTerminalMasses { get; }
        public double[] NTerminalMasses { get; }
        public double MonoisotopicMassIncludingFixedMods { get; }

        #endregion Public Properties

        #region Public Methods

        public override bool Equals(object obj)
        {
            var cp = obj as CompactPeptide;
            if (cp == null)
                return false;
            return (MonoisotopicMassIncludingFixedMods.Equals(cp.MonoisotopicMassIncludingFixedMods) && CTerminalMasses.SequenceEqual(cp.CTerminalMasses) &&
                NTerminalMasses.SequenceEqual(cp.NTerminalMasses));
        }

        public override int GetHashCode()
        {
            unchecked
            {
                var result = 0;
                foreach (double b in CTerminalMasses)
                    result = (result * 31) ^ b.GetHashCode();
                return result;
            }
        }

        public double[] ProductMassesMightHaveDuplicatesAndNaNs(List<ProductType> productTypes)
        {
            int massLen = 0;
            bool containsAdot = productTypes.Contains(ProductType.Adot);
            bool containsB = productTypes.Contains(ProductType.B);
            bool containsC = productTypes.Contains(ProductType.C);
            bool containsX = productTypes.Contains(ProductType.X);
            bool containsY = productTypes.Contains(ProductType.Y);
            bool containsZdot = productTypes.Contains(ProductType.Zdot);

            if (containsAdot)
                throw new NotImplementedException();
            if (containsB)
                massLen += NTerminalMasses.Length - 1;
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
            for (int j = 0; j < NTerminalMasses.Length; j++)
            {
                var hm = NTerminalMasses[j];
                if (containsB)
                {
                    if (j > 0)
                    {
                        massesToReturn[i] = hm;
                        i++;
                    }
                }
                if (containsC)
                {
                    massesToReturn[i] = hm + nitrogenAtomMonoisotopicMass + 3 * hydrogenAtomMonoisotopicMass;
                    i++;
                }
            }
            for (int j = 0; j < CTerminalMasses.Length; j++)
            {
                var hm = CTerminalMasses[j];
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

        #endregion Public Methods

        #region Private Methods

        private IEnumerable<double> ComputeFollowingFragmentMasses(PeptideWithSetModifications yyy, double prevMass, int oneBasedIndexToLookAt, int direction)
        {
            ModificationWithMass residue_variable_mod = null;
            do
            {
                prevMass += Residue.ResidueMonoisotopicMass[yyy[oneBasedIndexToLookAt - 1]];

                yyy.allModsOneIsNterminus.TryGetValue(oneBasedIndexToLookAt + 1, out residue_variable_mod);
                if (residue_variable_mod == null)
                {
                    yield return prevMass;
                }
                else if (residue_variable_mod.neutralLosses.Count == 1)
                {
                    prevMass += residue_variable_mod.monoisotopicMass - residue_variable_mod.neutralLosses.First();
                    yield return prevMass;
                }
                else
                {
                    foreach (double nl in residue_variable_mod.neutralLosses)
                    {
                        var theMass = prevMass + residue_variable_mod.monoisotopicMass - nl;
                        yield return theMass;
                        if ((direction == 1 && oneBasedIndexToLookAt + direction < yyy.Length) ||
                            (direction == -1 && oneBasedIndexToLookAt + direction > 1))
                            foreach (var nextMass in ComputeFollowingFragmentMasses(yyy, theMass, oneBasedIndexToLookAt + direction, direction))
                                yield return nextMass;
                    }
                    break;
                }
                oneBasedIndexToLookAt += direction;
            } while ((oneBasedIndexToLookAt > 1 && direction == -1) || (oneBasedIndexToLookAt < yyy.Length && direction == 1));
        }

        #endregion Private Methods

    }
}