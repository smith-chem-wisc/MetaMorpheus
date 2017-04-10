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

        public PsmClassic(PeptideWithSetModifications ps, string fileName, double scanRetentionTime, double scanPrecursorIntensity, double scanPrecursorMass, int scanNumber, int scanPrecursorNumber, int scanPrecursorCharge, int scanExperimentalPeaks, double totalIonCurrent, double scanPrecursorMZ, double score, int notch) : base(fileName, scanRetentionTime, scanPrecursorIntensity, scanPrecursorMass, scanNumber, scanPrecursorNumber, scanPrecursorCharge, scanExperimentalPeaks, totalIonCurrent, scanPrecursorMZ, score, notch)
        {
            this.ps = ps;
        }

        #endregion Public Constructors

        #region Public Methods

        public override CompactPeptide GetCompactPeptide(List<ModificationWithMass> variableModifications, List<ModificationWithMass> localizeableModifications, List<ModificationWithMass> fixedModifications)
        {
            if (compactPeptide == null)
                compactPeptide = new CompactPeptide(ps, variableModifications, localizeableModifications, fixedModifications);
            return compactPeptide;
        }

        #endregion Public Methods

        #region Internal Methods

        internal static bool FirstIsPreferable(PsmClassic psm, PsmClassic current_best_psm, List<ModificationWithMass> variableMods)
        {
            // Existed! Need to compare with old match
            if (Math.Abs(psm.score - current_best_psm.score) < 1e-9)
            {
                // Score is same, need to see if accepts and if prefer the new one
                return FirstIsPreferableWithoutScore(psm.ps, current_best_psm.ps, psm.scanPrecursorMass, variableMods);
            }
            if (psm.score > current_best_psm.score)
            {
                return true;
            }
            return false;
        }

        #endregion Internal Methods

        #region Private Methods

        private static bool FirstIsPreferableWithoutScore(PeptideWithSetModifications first, PeptideWithSetModifications second, double pm, List<ModificationWithMass> variableMods)
        {
            // If matches nicely to zero or to 1.003, and the other one doesn't, just accept it
            if ((Math.Abs(first.MonoisotopicMass - pm) < tolInDaForPreferringHavingMods || Math.Abs(first.MonoisotopicMass - pm + 1.003) < tolInDaForPreferringHavingMods || Math.Abs(first.MonoisotopicMass - pm + 2.0055) < tolInDaForPreferringHavingMods)
                && (Math.Abs(second.MonoisotopicMass - pm) > tolInDaForPreferringHavingMods && Math.Abs(second.MonoisotopicMass - pm + 1.003) > tolInDaForPreferringHavingMods && Math.Abs(second.MonoisotopicMass - pm + 2.0055) > tolInDaForPreferringHavingMods))
                return true;
            if ((Math.Abs(second.MonoisotopicMass - pm) < tolInDaForPreferringHavingMods || Math.Abs(second.MonoisotopicMass - pm + 1.003) < tolInDaForPreferringHavingMods || Math.Abs(second.MonoisotopicMass - pm + 2.0055) < tolInDaForPreferringHavingMods)
                && (Math.Abs(first.MonoisotopicMass - pm) > tolInDaForPreferringHavingMods && Math.Abs(first.MonoisotopicMass - pm + 1.003) > tolInDaForPreferringHavingMods && Math.Abs(first.MonoisotopicMass - pm + 2.0055) > tolInDaForPreferringHavingMods))
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

            return false;
        }

        #endregion Private Methods

    }
}