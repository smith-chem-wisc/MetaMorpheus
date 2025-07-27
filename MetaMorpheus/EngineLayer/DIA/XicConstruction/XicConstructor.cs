using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public abstract class XicConstructor
    {
        public Tolerance PeakFindingTolerance { get; set; } 
        public int MaxMissedScansAllowed { get; set; } 
        public double MaxPeakHalfWidth { get; set; }
        public int MinNumberOfPeaks { get; set; }

        public XicConstructor(Tolerance peakFindingTolerance, int maxMissedScansAllowed, double maxPeakHalfWidth, int minNumberOfPeaks)
        {
            this.PeakFindingTolerance = peakFindingTolerance;
            MaxMissedScansAllowed = maxMissedScansAllowed;
            MaxPeakHalfWidth = maxPeakHalfWidth;
            MinNumberOfPeaks = minNumberOfPeaks;
        }

        public abstract List<ExtractedIonChromatogram> GetAllXics(MsDataScan[] scans);
    }
}
