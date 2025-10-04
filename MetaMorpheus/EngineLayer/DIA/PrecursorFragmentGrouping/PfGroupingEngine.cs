using MassSpectrometry;
using System.Collections.Generic;

namespace EngineLayer.DIA
{
    /// <summary>
    /// PfGroupingEngine defines the process of grouping precursor and fragment XICs into PrecursorFragmentsGroup objects. 
    /// It should have a method that returns all PrecursorFragmentsGroup objects that can be found in a given set of precursor and fragment XICs.
    public abstract class PfGroupingEngine
    {
        public int MaxThreadsForGrouping { get; set; }
        public int MinFragmentCountForPfGroup { get; set; }
        public abstract List<PrecursorFragmentsGroup> PrecursorFragmentGrouping(List<ExtractedIonChromatogram> precursors, IEnumerable<ExtractedIonChromatogram> fragments);

        public static void FilterPfPairsByRank(List<PrecursorFragmentsGroup> pfGroups, int? precursorRankThreshold, int? fragmentRankThreshold)
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
                    pfGroup.PFpairs.RemoveAll(pf => pf.FragmentRank.Value > fragmentRankThreshold);
                }
            }
            if (precursorRankThreshold.HasValue)
            {
                foreach (var pfGroup in pfGroups)
                {
                    pfGroup.PFpairs.RemoveAll(pf => pf.PrecursorRank.Value > precursorRankThreshold);
                }
            }
        }
    }
}
