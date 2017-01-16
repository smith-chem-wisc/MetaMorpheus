using OldInternalLogic;
using System;
using System.Collections.Generic;

namespace InternalLogicEngineLayer
{
    public class ClassicSpectrumMatch : ParentSpectrumMatch
    {
        #region Public Fields

        public PeptideWithSetModifications ps;

        #endregion Public Fields

        #region Public Constructors

        public ClassicSpectrumMatch(double Score, PeptideWithSetModifications ps, double scanPrecursorMass, double scanPrecursorMZ, int scanNumber, double scanRT, int scanPrecursorCharge, int scanExperimentalPeaks, double TotalIonCurrent, double scanPrecursorIntensity, int spectraFileIndex)
        {
            this.ps = ps;
            this.Score = Score;
            this.scanPrecursorMass = scanPrecursorMass;

            this.scanPrecursorMZ = scanPrecursorMZ;
            this.scanNumber = scanNumber;
            this.scanPrecursorCharge = scanPrecursorCharge;
            this.scanRT = scanRT;
            this.scanPrecursorIntensity = scanPrecursorIntensity;
            this.scanExperimentalPeaks = scanExperimentalPeaks;
            this.TotalIonCurrent = TotalIonCurrent;
            this.spectraFileIndex = spectraFileIndex;
        }

        #endregion Public Constructors

        #region Public Properties

        public double scanPrecursorMZ { get; private set; }
        public double scanRT { get; private set; }
        public double scanPrecursorIntensity { get; private set; }
        public int scanExperimentalPeaks { get; private set; }
        public double TotalIonCurrent { get; private set; }
        public int spectraFileIndex { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public override CompactPeptide GetCompactPeptide(List<MorpheusModification> variableModifications, List<MorpheusModification> localizeableModifications)
        {
            if (compactPeptide == null)
                compactPeptide = new CompactPeptide(ps, variableModifications, localizeableModifications);
            return compactPeptide;
        }

        #endregion Public Methods

        #region Internal Methods

        internal static bool FirstIsPreferable(ClassicSpectrumMatch psm, ClassicSpectrumMatch current_best_psm)
        {
            // Existed! Need to compare with old match
            if (Math.Abs(psm.Score - current_best_psm.Score) < 1e-9)
            {
                // Score is same, need to see if accepts and if prefer the new one
                return FirstIsPreferableWithoutScore(psm.ps, current_best_psm.ps, psm.scanPrecursorMass);
            }
            if (psm.Score > current_best_psm.Score)
            {
                return true;
            }
            return false;
        }

        #endregion Internal Methods

        #region Private Methods

        private static bool FirstIsPreferableWithoutScore(PeptideWithSetModifications first, PeptideWithSetModifications second, double pm)
        {
            if (Math.Abs(first.MonoisotopicMass - pm) < 0.5 && Math.Abs(second.MonoisotopicMass - pm) > 0.5)
                return true;
            if (Math.Abs(first.MonoisotopicMass - pm) > 0.5 && Math.Abs(second.MonoisotopicMass - pm) < 0.5)
                return false;

            if (first.numVariableMods < second.numVariableMods)
                return true;

            return false;
        }

        #endregion Private Methods
    }
}