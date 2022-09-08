using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.MbrAnalysis
{
    public class TruePositiveDistribution : SpectralScoreDistribution
    {
        public TruePositiveDistribution(string outputFolder) : base(outputFolder)
        {
            SetDistributionType("TruePositive");
        }
    }
}
