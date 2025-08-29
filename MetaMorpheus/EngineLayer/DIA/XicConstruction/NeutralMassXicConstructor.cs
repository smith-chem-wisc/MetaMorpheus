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
        public double MinMass { get; set; } 
        public int MinCharge { get; set; }

        public NeutralMassXicConstructor(Tolerance peakFindingTolerance, int maxMissedScansAllowed, double maxPeakHalfWidth, int minNumberOfPeaks, DeconvolutionParameters deconParameters, double minMass = 0, int minCharge = 1, XicSpline? xicSpline = null) : base(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth, minNumberOfPeaks, xicSpline)
        {
            DeconParameters = deconParameters;
            MinMass = minMass;
            MinCharge = minCharge;
        }

        public override List<ExtractedIonChromatogram> GetAllXics(MsDataScan[] scans, MzRange isolationRange = null)
        {
            var neutralMassIndexingEngine = new MassIndexingEngine();
            if (neutralMassIndexingEngine.IndexPeaks(scans, DeconParameters, isolationRange, MinMass, MinCharge))
            {
                return neutralMassIndexingEngine.GetAllXics(PeakFindingTolerance, MaxMissedScansAllowed, MaxPeakHalfWidth, MinNumberOfPeaks, out var matchedPeaks);
            }
            else
            {
                throw new MetaMorpheusException("XIC construction failed.");
            }
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(BaseToString());
            sb.AppendLine($"NeutralMassXicConstructorSettings: ");
            sb.AppendLine($"DeconvolutionParameters: {DeconParameters}");
            sb.AppendLine($"MinMass: {MinMass}");
            sb.AppendLine($"MinCharge: {MinCharge}");
            return base.ToString();
        }
    }
}

