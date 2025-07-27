using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using MzLibUtil;

namespace EngineLayer.DIA
{
    public class NeutralMassXicConstructor : XicConstructor
    {
        public DeconvolutionParameters DeconParameters { get; set; } 
        public MzRange IsolationRange { get; set; } 
        public double MinMass { get; set; } 
        public int MinCharge { get; set; }

        public NeutralMassXicConstructor(Tolerance peakFindingTolerance, int maxMissedScansAllowed, double maxPeakHalfWidth, int minNumberOfPeaks, DeconvolutionParameters deconParameters, double minMass = 0, int minCharge = 1, MzRange isolationRange = null, XicSpline? xicSpline = null)
    : base(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth, minNumberOfPeaks, xicSpline)
        {
            DeconParameters = deconParameters;
            IsolationRange = isolationRange;
            MinMass = minMass;
            MinCharge = minCharge;
        }

        public override List<ExtractedIonChromatogram> GetAllXics(MsDataScan[] scans)
        {
            var neutralMassIndexingEngine = new MassIndexingEngine();
            if (neutralMassIndexingEngine.IndexPeaks(scans, DeconParameters, IsolationRange, MinMass, MinCharge))
            {
                return neutralMassIndexingEngine.GetAllXics(PeakFindingTolerance, MaxMissedScansAllowed, MaxPeakHalfWidth, MinNumberOfPeaks);
            }
            else
            {
                throw new MetaMorpheusException("XIC construction failed.");
            }
        }
    }
}

