using Chemistry;
using MetaMorpheus;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System;

namespace IndexSearchAndAnalyze
{
    public class ModernSpectrumMatch : ParentSpectrumMatch
    {
        public double ScoreFromSearch { get; private set; }
        public int spectraFileIndex { get; private set; }
        public double scanRT { get; private set; }
        public double scanPrecursorMZ { get; private set; }
        public double scanPrecursorIntensity { get; private set; }
        public int scanExperimentalPeaks { get; private set; }
        public double TotalIonCurrent { get; private set; }

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
            this.ScoreFromSearch = score;
            this.spectraFileIndex = spectraFileIndex;
            this.compactPeptide = theBestPeptide;
        }

        internal override CompactPeptide GetCompactPeptide(List<MorpheusModification> variableModifications, List<MorpheusModification> localizeableModifications)
        {
            return compactPeptide;
        }
    }
}