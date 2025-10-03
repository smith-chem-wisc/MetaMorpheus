using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using MassSpectrometry;
using System.Collections.Concurrent;
using System.Linq;

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
        /// Given a list of precursor XICs and all eligible fragment XICs, loop over each precursor and group fragments with the precursor based on the grouping criteria.
        public override IEnumerable<PrecursorFragmentsGroup> PrecursorFragmentGrouping(List<ExtractedIonChromatogram> precursors, List<ExtractedIonChromatogram> fragments)
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
                        var pfGroup = GroupFragmentsForOnePrecursor(precursor, fragmentsInRange, OverlapThreshold, CorrelationThreshold, MinFragmentCountForPfGroup);
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
        /// </summary>
        public static PrecursorFragmentsGroup GroupFragmentsForOnePrecursor(ExtractedIonChromatogram precursorXic, List<ExtractedIonChromatogram> fragmentXics, double overlapThreshold, double correlationThreshold, int minFragmentCountForGrouping)
        {
            var pfPairs = new List<PrecursorFragmentPair>();
            foreach (var fragmentXic in fragmentXics)
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
            if (pfPairs.Count > minFragmentCountForGrouping)
            {
                var pfGroup = new PrecursorFragmentsGroup(precursorXic, pfPairs);
                return pfGroup;
            }
            return null;
        }

        protected static SortedDictionary<double, List<ExtractedIonChromatogram>> BuildApexSortedXics(List<ExtractedIonChromatogram> xics)
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
    }
}
