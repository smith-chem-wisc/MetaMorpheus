using Chemistry;
using OldInternalLogic;
using System.Collections.Generic;

namespace InternalLogicEngineLayer
{
    public class ModernSpectrumMatch : ParentSpectrumMatch
    {

        #region Public Constructors

        public ModernSpectrumMatch(double scanPrecursorMZ, int scanNumber, double scanRT, int scanPrecursorCharge, int scanExperimentalPeaksCount, double totalIonCurrent, double precursorIntensity, int spectraFileIndex, CompactPeptide theBestPeptide, double score)
        {
            this.ScanPrecursorMZ = scanPrecursorMZ;
            this.scanNumber = scanNumber;
            this.scanPrecursorCharge = scanPrecursorCharge;
            this.ScanRT = scanRT;
            scanPrecursorMass = scanPrecursorMZ.ToMass(scanPrecursorCharge);
            ScanPrecursorIntensity = precursorIntensity;
            ScanExperimentalPeaks = scanExperimentalPeaksCount;
            TotalIonCurrent = totalIonCurrent;
            Score = score;
            this.SpectraFileIndex = spectraFileIndex;
            compactPeptide = theBestPeptide;
        }

        #endregion Public Constructors

        #region Public Properties

        public int SpectraFileIndex { get; private set; }
        public double ScanRT { get; private set; }
        public double ScanPrecursorMZ { get; private set; }
        public double ScanPrecursorIntensity { get; private set; }
        public int ScanExperimentalPeaks { get; private set; }
        public double TotalIonCurrent { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public override CompactPeptide GetCompactPeptide(List<MorpheusModification> variableModifications, List<MorpheusModification> localizeableModifications)
        {
            return compactPeptide;
        }

        #endregion Public Methods

    }
}