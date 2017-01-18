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

        public ClassicSpectrumMatch(double score, PeptideWithSetModifications ps, double scanPrecursorMass, double scanPrecursorMZ, int scanNumber, double scanRT, int scanPrecursorCharge, int scanExperimentalPeaks, double totalIonCurrent, double scanPrecursorIntensity, int spectraFileIndex)
        {
            this.ps = ps;
            this.Score = score;
            this.scanPrecursorMass = scanPrecursorMass;

            this.ScanPrecursorMZ = scanPrecursorMZ;
            this.scanNumber = scanNumber;
            this.scanPrecursorCharge = scanPrecursorCharge;
            this.ScanRT = scanRT;
            this.ScanPrecursorIntensity = scanPrecursorIntensity;
            this.ScanExperimentalPeaks = scanExperimentalPeaks;
            this.TotalIonCurrent = totalIonCurrent;
            this.SpectraFileIndex = spectraFileIndex;
        }

        #endregion Public Constructors

        #region Public Properties

        public double ScanPrecursorMZ { get; private set; }
        public double ScanRT { get; private set; }
        public double ScanPrecursorIntensity { get; private set; }
        public int ScanExperimentalPeaks { get; private set; }
        public double TotalIonCurrent { get; private set; }
        public int SpectraFileIndex { get; private set; }

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

            if (first.NumVariableMods < second.NumVariableMods)
                return true;

            return false;
        }

        #endregion Private Methods
    }
}