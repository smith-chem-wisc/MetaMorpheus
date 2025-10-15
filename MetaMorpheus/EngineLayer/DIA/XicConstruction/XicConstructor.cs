using MassSpectrometry;
using MzLibUtil;
using System.Collections.Generic;

namespace EngineLayer.DIA
{
    /// <summary>
    /// XicConstructor defines the process of constructing extracted ion chromatograms (XICs) for a given set of scans. It could operate on spectral peaks
    /// or deconvoluted masses depending on the implementation. Any XicConstructor should have a method that returns all XICs that can be found in the set of scans.
    /// </summary>
    public abstract class XicConstructor
    {
        public Tolerance PeakFindingTolerance { get; set; } 
        public int MaxMissedScansAllowed { get; set; } 
        public double MaxPeakHalfWidth { get; set; }
        public int MinNumberOfPeaks { get; set; }
        public XicSpline? XicSpline { get; set; }

        public XicConstructor(Tolerance peakFindingTolerance, int maxMissedScansAllowed, double maxPeakHalfWidth, int minNumberOfPeaks, XicSpline? xicSpline = null)
        {
            this.PeakFindingTolerance = peakFindingTolerance;
            MaxMissedScansAllowed = maxMissedScansAllowed;
            MaxPeakHalfWidth = maxPeakHalfWidth;
            MinNumberOfPeaks = minNumberOfPeaks;
            XicSpline = xicSpline;
        }

        public abstract List<ExtractedIonChromatogram> GetAllXics(MsDataScan[] scans, MzRange isolationRange = null);
    }
}
