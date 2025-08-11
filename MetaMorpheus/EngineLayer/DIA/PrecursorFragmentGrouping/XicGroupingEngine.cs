using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using MassSpectrometry;
using System.Collections.Concurrent;

namespace EngineLayer.DIA
{
    /// <summary>
    /// XicGroupingEngine is a specific implementation of PfGroupingEngine that groups precursor and fragment XICs into PrecursorFragmentsGroup objects 
    /// based on specified criteria including apex RT tolerance, overlap threshold, and correlation threshold.
    /// </summary>
    public class XicGroupingEngine : PfGroupingEngine
    {
        public float ApexRTTolerance { get; set; }
        public double OverlapThreshold { get; set; }
        public double CorrelationThreshold { get; set; }

        public XicGroupingEngine(float apexRTTolerance, double overlapThreshold, double correlationThreshold, int maxThreadsForGrouping = 1, int minFragmentCountForGrouping = 0)
        {
            ApexRTTolerance = apexRTTolerance;
            OverlapThreshold = overlapThreshold;
            CorrelationThreshold = correlationThreshold;
            MaxThreadsForGrouping = maxThreadsForGrouping;
            MinFragmentCountForPfGroup = minFragmentCountForGrouping;
        }

        /// <summary>
        /// Given a list of precursor XICs and all eligibile fragment XICs, loop over each precursor and group fragments with the precursor based on the grouping criteria.
        /// <summary>
        public override List<PrecursorFragmentsGroup> PrecursorFragmentGrouping(List<ExtractedIonChromatogram> precursors, List<ExtractedIonChromatogram> fragments)
        {
            var pfGroups = new List<PrecursorFragmentsGroup>();

            Parallel.ForEach(Partitioner.Create(0, precursors.Count), new ParallelOptions { MaxDegreeOfParallelism = MaxThreadsForGrouping },
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        var precursor = precursors[i];
                        var pfGroup = GroupFragmentsForOnePrecursor(precursor, fragments, ApexRTTolerance, OverlapThreshold, CorrelationThreshold, MinFragmentCountForPfGroup);
                        if (pfGroup != null)
                        {
                            lock (pfGroups)
                                pfGroups.Add(pfGroup);
                        }
                    }
                });
            return pfGroups;
        }

        /// <summary>
        /// Given one precursor XIC and all eligibile fragment XICs, select fragments that meet the grouping criteria and construct a precursor-fragment group for this precursor.
        /// <summary>
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
                        double correlation = PrecursorFragmentsGroup.CalculateXicCorrelation(precursorXic, fragmentXic);
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
