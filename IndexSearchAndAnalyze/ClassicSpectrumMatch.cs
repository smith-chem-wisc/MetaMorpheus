using System;
using System.Collections.Generic;
using MassSpectrometry;
using Spectra;
using System.Linq;

namespace IndexSearchAndAnalyze
{
    internal class ClassicSpectrumMatch
    {
        private double[] hehe;
        public CompactPeptide ps;

        public ClassicSpectrumMatch(double score, CompactPeptide ps, double[] hehe, double precursorMass, double scanPrecursorMZ, int scanNumber, double scanRT, int scanPrecursorCharge, int scanExperimentalPeaksCount, double totalIonCurrent, double precursorIntensity, int spectraFileIndex)
        {
            this.ps = ps;
            this.hehe = hehe;
            this.score = score;
            this.precursorMass = precursorMass;

            this.scanPrecursorMZ = scanPrecursorMZ;
            this.scanNumber = scanNumber;
            this.scanPrecursorCharge = scanPrecursorCharge;
            this.scanRT = scanRT;
            this.scanPrecursorIntensity = precursorIntensity;
            this.scanExperimentalPeaks = scanExperimentalPeaksCount;
            this.TotalIonCurrent = totalIonCurrent;
            this.spectraFileIndex = spectraFileIndex;
        }

        public double precursorMass { get; private set; }
        public double score { get; private set; }
        public double scanPrecursorMZ { get; private set; }
        public int scanNumber { get; private set; }
        public int scanPrecursorCharge { get; private set; }
        public double scanRT { get; private set; }
        public double scanPrecursorIntensity { get; private set; }
        public int scanExperimentalPeaks { get; private set; }
        public double TotalIonCurrent { get; private set; }
        public double ScoreFromSearch { get; private set; }
        public int spectraFileIndex { get; private set; }

        internal static bool FirstIsPreferable(ClassicSpectrumMatch psm, ClassicSpectrumMatch current_best_psm)
        {
            // Existed! Need to compare with old match
            if (Math.Abs(psm.score - current_best_psm.score) < 1e-9)
            {
                // Score is same, need to see if accepts and if prefer the new one
                return ModernSearchEngine.FirstIsPreferableWithoutScore(psm.ps, current_best_psm.ps, psm.precursorMass);
            }
            else if (psm.score > current_best_psm.score)
            {
                return true;
            }
            return false;
        }
    }
}