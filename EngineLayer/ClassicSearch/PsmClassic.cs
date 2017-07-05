using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.ClassicSearch
{
    public class PsmClassic : PsmParent
    {
        #region Public Fields

        public PeptideWithSetModifications ps;

        #endregion Public Fields

        #region Private Fields

        private const double tolInDaForPreferringHavingMods = 0.03;
        private CompactPeptide compactPeptide;

        #endregion Private Fields

        #region Public Constructors

        public PsmClassic(PeptideWithSetModifications ps, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan) : base(notch, score, scanIndex, scan, ps.MonoisotopicMass)
        {
            this.ps = ps;
        }

        #endregion Public Constructors

        #region Public Methods

        public override CompactPeptide GetCompactPeptide(Dictionary<ModificationWithMass, ushort> modsDictionary)
        {
            if (compactPeptide == null)
                compactPeptide = new CompactPeptide(ps, modsDictionary);
            return compactPeptide;
        }

        #endregion Public Methods

        #region Internal Methods

        internal static bool? FirstIsPreferable(PsmClassic firstPsm, PsmClassic secondPsm, List<ModificationWithMass> variableMods)
        {
            // Existed! Need to compare with old match
            if (Math.Abs(firstPsm.Score - secondPsm.Score) < 1e-9)
            {
                // Score is same, need to see if accepts and if prefer the new one
                var first = firstPsm.ps;
                var second = secondPsm.ps;
                var pm = firstPsm.ScanPrecursorMass;

                // Prefer to be at zero rather than fewer modifications
                if ((Math.Abs(first.MonoisotopicMass - pm) < tolInDaForPreferringHavingMods)
                    && (Math.Abs(second.MonoisotopicMass - pm) > tolInDaForPreferringHavingMods))
                    return true;
                if ((Math.Abs(second.MonoisotopicMass - pm) < tolInDaForPreferringHavingMods)
                    && (Math.Abs(first.MonoisotopicMass - pm) > tolInDaForPreferringHavingMods))
                    return false;

                int firstVarMods = first.allModsOneIsNterminus.Count(b => variableMods.Contains(b.Value));
                int secondVarMods = second.allModsOneIsNterminus.Count(b => variableMods.Contains(b.Value));

                // Want the lowest number of localizeable mods!!! Even at the expense of many variable and fixed mods.
                if (first.NumMods - firstVarMods - first.numFixedMods < second.NumMods - secondVarMods - second.numFixedMods)
                    return true;
                if (first.NumMods - firstVarMods - first.numFixedMods > second.NumMods - secondVarMods - second.numFixedMods)
                    return false;

                // If have same number of localizeable mods, pick the lowest number of localizeable + variable mods
                if (first.NumMods - first.numFixedMods < second.NumMods - second.numFixedMods)
                    return true;
                if (first.NumMods - first.numFixedMods > second.NumMods - second.numFixedMods)
                    return false;

                // If have same numbers of localizeable and variable mods, prefer not to have substitutions and removals!
                int firstNumRazor = first.allModsOneIsNterminus.Count(b => b.Value.modificationType.Equals("substitution") || b.Value.modificationType.Equals("missing") || b.Value.modificationType.Equals("trickySubstitution"));
                int secondNumRazor = second.allModsOneIsNterminus.Count(b => b.Value.modificationType.Equals("substitution") || b.Value.modificationType.Equals("missing") || b.Value.modificationType.Equals("trickySubstitution"));

                if (firstNumRazor < secondNumRazor)
                    return true;
                if (firstNumRazor > secondNumRazor)
                    return false;

                return null;
            }
            if (firstPsm.Score > secondPsm.Score)
            {
                return true;
            }
            return false;
        }

        #endregion Internal Methods
    }
}