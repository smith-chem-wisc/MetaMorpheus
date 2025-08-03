using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using System.Collections.Concurrent;
using ThermoFisher.CommonCore.Data.Business;
using System.Text.RegularExpressions;

namespace EngineLayer.DIA
{
    /// <summary>
    /// XicGroupingEngine is a specific implementation of PfGroupingEngine that groups precursor and fragment XICs into PrecursorFragmentsGroup objects 
    /// based on specified criteria including apex RT tolerance, overlap threshold, and correlation threshold.
    /// <summary>
    public class XicGroupingEngine : PfGroupingEngine
    {
        public float ApexRTTolerance { get; set; }
        public double OverlapThreshold { get; set; }
        public double CorrelationThreshold { get; set; }
        public int PrecursorRankThreshold { get; set; } 
        public int FragmentRankThreshold { get; set; }

        public XicGroupingEngine(float apexRTTolerance, double overlapThreshold, double correlationThreshold, int maxThreadsForGrouping = 1, int minFragmentCountForGrouping = 0, int precursorRankThreshold = 1000, int fragmentRankThreshold = 1000)
        {
            ApexRTTolerance = apexRTTolerance;
            OverlapThreshold = overlapThreshold;
            CorrelationThreshold = correlationThreshold;
            MaxThreadsForGrouping = maxThreadsForGrouping;
            MinFragmentCountForPfGroup = minFragmentCountForGrouping;
            PrecursorRankThreshold = precursorRankThreshold;
            FragmentRankThreshold = fragmentRankThreshold;
        }

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

            //Filter precursor-fragment pairs by rank
            FilterPfPairsByRank(pfGroups, PrecursorRankThreshold, FragmentRankThreshold);
            //Remove groups with insufficient fragment pairs
            pfGroups.RemoveAll(g => g.PFpairs.Count < MinFragmentCountForPfGroup);

            return pfGroups;
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

        public static void FilterPfPairsByRank(List<PrecursorFragmentsGroup> pfGroups, int? fragmentRankThreshold, int? precursorRankThreshold)
        {
            //Rank precursors for all precursor-fragment pairs in all precursor-fragment groups
            if (precursorRankThreshold.HasValue)
            {
                PrecursorFragmentPair.SetPrecursorRankForPfPairs(pfGroups.SelectMany(g => g.PFpairs));
            }
            //Rank fragments for all precursor-fragment pairs within each precursor-fragment group
            if (fragmentRankThreshold.HasValue)
            {
                foreach (var pfGroup in pfGroups)
                {
                    pfGroup.SetFragmentRankForPfPairs();
                    pfGroup.PFpairs.RemoveAll(pf => pf.FragmentRank.Value < fragmentRankThreshold);
                }
            }
            if (precursorRankThreshold.HasValue)
            {
                foreach (var pfGroup in pfGroups)
                {
                    pfGroup.PFpairs.RemoveAll(pf => pf.PrecursorRank.Value < precursorRankThreshold);
                }
            }
        }
    }
}
