using MassSpectrometry;
using MzLibUtil;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ThermoFisher.CommonCore.Data.Business;

namespace EngineLayer.DIA
{
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
