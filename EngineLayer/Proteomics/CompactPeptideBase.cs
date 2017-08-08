using Chemistry;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    [Serializable]
    public abstract class CompactPeptideBase
    {
        #region Protected Fields

        protected static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        protected static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        protected static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;
        protected static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        #endregion Protected Fields

        #region Public Properties

        public double[] CTerminalMasses { get; protected set; }
        public double[] NTerminalMasses { get; protected set; }
        public double MonoisotopicMassIncludingFixedMods { get; protected set; }

        #endregion Public Properties

        #region Public Methods

        public override bool Equals(object obj)
        {
            var cp = obj as CompactPeptideBase;
            if (cp == null)
                return false;
            if (CTerminalMasses.Count() == 1 || cp.CTerminalMasses.Count() == 1)
            {
                return (
                    ((double.IsNaN(MonoisotopicMassIncludingFixedMods) && double.IsNaN(cp.MonoisotopicMassIncludingFixedMods)) || Math.Abs(MonoisotopicMassIncludingFixedMods - cp.MonoisotopicMassIncludingFixedMods) < 1e-7)
                    && NTerminalMasses.SequenceEqual(cp.NTerminalMasses)
                    );
            }
            else if (NTerminalMasses.Count() == 1 || cp.NTerminalMasses.Count() == 1)
            {
                return (
                    ((double.IsNaN(MonoisotopicMassIncludingFixedMods) && double.IsNaN(cp.MonoisotopicMassIncludingFixedMods)) || Math.Abs(MonoisotopicMassIncludingFixedMods - cp.MonoisotopicMassIncludingFixedMods) < 1e-7)
                    && CTerminalMasses.SequenceEqual(cp.CTerminalMasses)
                    );
            }
            else
            {
                return (
                    ((double.IsNaN(MonoisotopicMassIncludingFixedMods) && double.IsNaN(cp.MonoisotopicMassIncludingFixedMods)) || Math.Abs(MonoisotopicMassIncludingFixedMods - cp.MonoisotopicMassIncludingFixedMods) < 1e-7)
                    && CTerminalMasses.SequenceEqual(cp.CTerminalMasses)
                    && NTerminalMasses.SequenceEqual(cp.NTerminalMasses)
                    );
            }
        }

        public override int GetHashCode()
        {
            unchecked
            {
                var result = 0;
                if (CTerminalMasses.Count() == 1)
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

        #endregion Public Methods

        #region Protected Methods

        protected IEnumerable<double> ComputeFollowingFragmentMasses(PeptideWithSetModifications yyy, double prevMass, int oneBasedIndexToLookAt, int direction)
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

        #endregion Protected Methods
    }
}