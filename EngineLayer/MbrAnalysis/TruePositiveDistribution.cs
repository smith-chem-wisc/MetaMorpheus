using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.MbrAnalysis
{
    public class TruePositiveDistribution : SpectralScoreDistribution
    {
        public List<SpectralComparison> SpectralComparisonList { get; private set; }
        public string OutputFolder { get; set; }
        public string DistributionType { get; private set; }

        public TruePositiveDistribution(string outputFolder) : base(outputFolder)
        {
            DistributionType = "TruePositive";
        }
    }
}
