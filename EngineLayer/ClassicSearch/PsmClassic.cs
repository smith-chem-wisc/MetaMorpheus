using Proteomics;
using System;
using System.Collections.Generic;

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

        public PsmClassic(PeptideWithSetModifications ps, string fileName, double scanRetentionTime, double scanPrecursorIntensity, double scanPrecursorMass, int scanNumber, int scanPrecursorCharge, int scanExperimentalPeaks, double totalIonCurrent, double scanPrecursorMZ, double score, int notch) : base(fileName, scanRetentionTime, scanPrecursorIntensity, scanPrecursorMass, scanNumber, scanPrecursorCharge, scanExperimentalPeaks, totalIonCurrent, scanPrecursorMZ, score, notch)
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

        internal static bool FirstIsPreferable(PsmClassic psm, PsmClassic current_best_psm)
        {
            // Existed! Need to compare with old match
            if (Math.Abs(psm.score - current_best_psm.score) < 1e-9)
            {
                // Score is same, need to see if accepts and if prefer the new one
                return FirstIsPreferableWithoutScore(psm.ps, current_best_psm.ps, psm.scanPrecursorMass);
            }
            if (psm.score > current_best_psm.score)
            {
                return true;
            }
            return false;
        }

        #endregion Internal Methods

        #region Private Methods

        private static bool FirstIsPreferableWithoutScore(PeptideWithSetModifications first, PeptideWithSetModifications second, double pm)
        {
            if (Math.Abs(first.MonoisotopicMass - pm) < tolInDaForPreferringHavingMods && Math.Abs(second.MonoisotopicMass - pm) > tolInDaForPreferringHavingMods)
                return true;
            if (Math.Abs(first.MonoisotopicMass - pm) > tolInDaForPreferringHavingMods && Math.Abs(second.MonoisotopicMass - pm) < tolInDaForPreferringHavingMods)
                return false;
            if (first.NumMods < second.NumMods)
                return true;

            return false;
        }

        #endregion Private Methods

    }
}