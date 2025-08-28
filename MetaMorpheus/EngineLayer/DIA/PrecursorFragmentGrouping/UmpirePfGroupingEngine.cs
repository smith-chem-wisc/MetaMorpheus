using MassSpectrometry;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Configuration;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class UmpirePfGroupingEngine : XicGroupingEngine
    {
        public int NoPointPerInterval { get; set; }
        public UmpirePfGroupingEngine(int noPointPerInterval, float apexRTTolerance, double overlapThreshold, double correlationThreshold, int maxThreadsForGrouping = 1, int minFragmentCountForGrouping = 0, int? precursorRankThreshold = null, int? fragmentRankThreshold = null)
            : base(apexRTTolerance, overlapThreshold, correlationThreshold, maxThreadsForGrouping, minFragmentCountForGrouping, precursorRankThreshold, fragmentRankThreshold)
        {
            NoPointPerInterval = noPointPerInterval;
        }

        public override PrecursorFragmentsGroup GroupFragmentsForOnePrecursor(ExtractedIonChromatogram precursorXic, IEnumerable<ExtractedIonChromatogram> fragmentXics)
        {
            var pfPairs = new List<PrecursorFragmentPair>();
            foreach (var fragmentXic in fragmentXics)
            {
                double overlap = PrecursorFragmentsGroup.CalculateXicOverlapRatio_Umpire(precursorXic, fragmentXic);
                if (overlap >= OverlapThreshold)
                {
                    double correlation = PrecursorFragmentsGroup.CalculateXicCorrXYData_Umpire(precursorXic, fragmentXic, NoPointPerInterval);
                    if (correlation >= CorrelationThreshold)
                    {
                        var pfPair = new PrecursorFragmentPair(precursorXic, fragmentXic, correlation, overlap);
                        pfPairs.Add(pfPair);
                    }
                }
            }
            if (pfPairs.Count > MinFragmentCountForPfGroup)
            {
                var pfGroup = new PrecursorFragmentsGroup(precursorXic, pfPairs);
                return pfGroup;
            }
            return null;
        }
    }
}
