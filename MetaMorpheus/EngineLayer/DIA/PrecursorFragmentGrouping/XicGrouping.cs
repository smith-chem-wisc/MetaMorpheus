using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Security.Policy;

namespace EngineLayer.DIA
{
    public class XicGrouping : PfGroupingEngine
    {
        public float ApexRTTolerance { get; set; }
        public double OverlapThreshold { get; set; }
        public double CorrelationThreshold { get; set; }
        public int MaxThreadsForGrouping { get; set; } 
        public int MinFragmentCountForPfGroup { get; set; } 

        public XicGrouping(float apexRTTolerance, double overlapThreshold, double correlationThreshold, int maxThreadsForGrouping = 1, int minFragmentCountForGrouping = 0)
        {
            ApexRTTolerance = apexRTTolerance;
            OverlapThreshold = overlapThreshold;
            CorrelationThreshold = correlationThreshold;
            MaxThreadsForGrouping = maxThreadsForGrouping;
            MinFragmentCountForPfGroup = minFragmentCountForGrouping;
        }

        public override List<PrecursorFragmentsGroup> PrecursorFragmentGrouping(List<ExtractedIonChromatogram> precursors, List<ExtractedIonChromatogram> fragments)
        {
            var precursorGroups = new List<PrecursorFragmentsGroup>();
            foreach (var precursorXic in precursors)
            {
                var pfGroup = GroupFragmentsForOnePrecursor(precursorXic, fragments, ApexRTTolerance, OverlapThreshold, CorrelationThreshold, MinFragmentCountForPfGroup);
                if (pfGroup != null)
                {
                    precursorGroups.Add(pfGroup);
                }
            }
            return precursorGroups;
        }

        public static PrecursorFragmentsGroup GroupFragmentsForOnePrecursor(ExtractedIonChromatogram precursorXic, List<ExtractedIonChromatogram> fragmentXics, float apexRtTolerance, double overlapThreshold, double correlationThreshold, int minFragmentCountForGrouping)
        {
            var pfPairs = new List<PrecursorFragmentPair>();
            foreach (var fragmentXic in fragmentXics)
            {
                if (Math.Abs(fragmentXic.ApexRT - precursorXic.ApexRT) <= apexRtTolerance)
                {
                    double overlap = PrecursorFragmentsGroup.CalculateXicOverlapRatio(precursorXic, fragmentXic);
                    if (overlap >= overlapThreshold)
                    {
                        double correlation = PrecursorFragmentsGroup.CalculateXicCorrelationXYData(precursorXic, fragmentXic);
                        if (correlation >= correlationThreshold)
                        {
                            var pfPair = new PrecursorFragmentPair(precursorXic, fragmentXic, correlation, overlap);
                            pfPairs.Add(pfPair);
                        }
                    }
                }
            }
            if (pfPairs.Count > minFragmentCountForGrouping)
            {
                var pfGroup = new PrecursorFragmentsGroup(precursorXic, pfPairs);
                return pfGroup;
            }
            return null;
        }
    }
}
