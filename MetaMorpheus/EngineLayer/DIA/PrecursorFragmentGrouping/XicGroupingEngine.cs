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
        public int? PrecursorRankThreshold { get; set; } 
        public int? FragmentRankThreshold { get; set; }

        public XicGroupingEngine(float apexRTTolerance, double overlapThreshold, double correlationThreshold, int maxThreadsForGrouping = 1, int minFragmentCountForGrouping = 0, int? precursorRankThreshold = null, int? fragmentRankThreshold = null)
        {
            ApexRTTolerance = apexRTTolerance;
            OverlapThreshold = overlapThreshold;
            CorrelationThreshold = correlationThreshold;
            MaxThreadsForGrouping = maxThreadsForGrouping;
            MinFragmentCountForPfGroup = minFragmentCountForGrouping;
            PrecursorRankThreshold = precursorRankThreshold;
            FragmentRankThreshold = fragmentRankThreshold;
        }

        public override List<PrecursorFragmentsGroup> PrecursorFragmentGrouping(List<ExtractedIonChromatogram> precursors, IEnumerable<ExtractedIonChromatogram> fragments)
        {
            var pfGroups = new List<PrecursorFragmentsGroup>();
            var apexSortedFragmentXics = BuildApexSortedXics(fragments);

            Parallel.ForEach(Partitioner.Create(0, precursors.Count), new ParallelOptions { MaxDegreeOfParallelism = MaxThreadsForGrouping },
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        var precursor = precursors[i];
                        var fragmentsInRange = GetXicsInRange(apexSortedFragmentXics, precursor.ApexRT, ApexRTTolerance);
                        var pfGroup = GroupFragmentsForOnePrecursor(precursor, fragmentsInRange);
                        if (pfGroup != null)
                        {
                            lock (pfGroups)
                                pfGroups.Add(pfGroup);
                        }
                    }
                });

            //Filter precursor-fragment pairs by rank if rank thresholds are set
            FilterPfPairsByRank(pfGroups, PrecursorRankThreshold, FragmentRankThreshold);

            //Remove groups with insufficient fragment pairs after filtering
            pfGroups.RemoveAll(g => g.PFpairs.Count < MinFragmentCountForPfGroup);

            return pfGroups;
        }

        public virtual PrecursorFragmentsGroup GroupFragmentsForOnePrecursor(ExtractedIonChromatogram precursorXic, IEnumerable<ExtractedIonChromatogram> fragmentXics)
        {
            return GroupFragmentsForOnePrecursor(precursorXic, fragmentXics, ApexRTTolerance, OverlapThreshold, CorrelationThreshold, MinFragmentCountForPfGroup);
        }

        public static PrecursorFragmentsGroup GroupFragmentsForOnePrecursor(ExtractedIonChromatogram precursorXic, IEnumerable<ExtractedIonChromatogram> fragmentXics, float apexRtTolerance, double overlapThreshold, double correlationThreshold, int minFragmentCountForGrouping)
        {
            var pfPairs = new List<PrecursorFragmentPair>();
            foreach (var fragmentXic in fragmentXics)
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
            if (pfPairs.Count > minFragmentCountForGrouping)
            {
                var pfGroup = new PrecursorFragmentsGroup(precursorXic, pfPairs);
                return pfGroup;
            }
            return null;
        }

        protected static SortedDictionary<double, List<ExtractedIonChromatogram>> BuildApexSortedXics(IEnumerable<ExtractedIonChromatogram> xics)
        {
            var tree = new SortedDictionary<double, List<ExtractedIonChromatogram>>();
            foreach (var xic in xics)
            {
                double roundedApexRt = Math.Round(xic.ApexRT, 2);
                if (!tree.ContainsKey(roundedApexRt))
                    tree[roundedApexRt] = new List<ExtractedIonChromatogram>();
                tree[roundedApexRt].Add(xic);
            }
            return tree;
        }

        protected static List<ExtractedIonChromatogram> GetXicsInRange(SortedDictionary<double, List<ExtractedIonChromatogram>> apexSortedFragmentXics, double targetApexRt, double rtTolerance)
        {
            double minRt = Math.Round(targetApexRt - rtTolerance, 2);
            double maxRt = Math.Round(targetApexRt + rtTolerance, 2);

            var subset = apexSortedFragmentXics.Where(kvp => kvp.Key >= minRt && kvp.Key <= maxRt).SelectMany(kvp => kvp.Value).ToList();
            return subset;
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine($"XicGroupingEngineSettings: ");
            sb.AppendLine($"ApexRTTolerance: {ApexRTTolerance}");
            sb.AppendLine($"OverlapThreshold: {OverlapThreshold}");
            sb.AppendLine($"CorrelationThreshold: {CorrelationThreshold}");
            if (PrecursorRankThreshold.HasValue)
            {
                sb.AppendLine($"PrecursorRankThreshold: {PrecursorRankThreshold}");
            }
            else
            {
                sb.AppendLine($"PrecursorRankThreshold: None");
            }
            if (FragmentRankThreshold.HasValue)
            {
                sb.AppendLine($"FragmentRankThreshold: {FragmentRankThreshold}");
            }
            else
            {
                sb.AppendLine($"FragmentRankThreshold: None");
            }
            return sb.ToString();
        }
    }
}
