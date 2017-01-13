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
            this.scanPrecursorMZ = scanPrecursorMZ;
            this.scanNumber = scanNumber;
            this.scanPrecursorCharge = scanPrecursorCharge;
            this.scanRT = scanRT;
            this.scanPrecursorMass = scanPrecursorMZ.ToMass(scanPrecursorCharge);
            this.scanPrecursorIntensity = precursorIntensity;
            this.scanExperimentalPeaks = scanExperimentalPeaksCount;
            this.TotalIonCurrent = totalIonCurrent;
            this.Score = score;
            this.spectraFileIndex = spectraFileIndex;
            this.compactPeptide = theBestPeptide;
        }

        #endregion Public Constructors

        #region Public Properties

        public double ScoreFromSearch { get; private set; }
        public int spectraFileIndex { get; private set; }
        public double scanRT { get; private set; }
        public double scanPrecursorMZ { get; private set; }
        public double scanPrecursorIntensity { get; private set; }
        public int scanExperimentalPeaks { get; private set; }
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