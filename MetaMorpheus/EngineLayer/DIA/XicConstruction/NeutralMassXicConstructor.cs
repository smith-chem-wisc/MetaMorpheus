using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using MzLibUtil;

namespace EngineLayer.DIA.XicConstruction
{
    public class NeutralMassXicConstructor : XicConstructor
    {
        public DeconvolutionParameters DeconParameters { get; set; } 
        public double MinMass { get; set; } 
        public double MaxMass { get; set; }
        public int MinCharge { get; set; }
        public int MaxCharge { get; set; }
        
        public NeutralMassXicConstructor(Tolerance peakFindingTolerance, int maxMissedScansAllowed, double maxPeakHalfWidth, int minNumberOfPeaks, XicSpline xicSpline)
    : base(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth, minNumberOfPeaks, xicSpline)
        {
        }

        public override List<ExtractedIonChromatogram> GetAllXics(MsDataScan[] scans)
        {
            var neutralMassIndexingEngine = MassIndexingEngine.InitializeMassIndexingEngine(scans, DeconParameters);
            return neutralMassIndexingEngine.GetAllXics(PeakFindingTolerance, MaxMissedScansAllowed, MaxPeakHalfWidth, MinNumberOfPeaks);
        }
    }
}

