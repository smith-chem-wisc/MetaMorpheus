using MassSpectrometry;
using System.Collections.Generic;

namespace EngineLayer.DIA
{
    /// <summary>
    /// PfGroupingEngine defines the process of grouping precursor and fragment XICs into PrecursorFragmentsGroup objects. 
    /// It should have a method that returns all PrecursorFragmentsGroup objects that can be found in a given set of precursor and fragment XICs.
    /// <summary>
    public abstract class PfGroupingEngine
    {
        public int MaxThreadsForGrouping { get; set; }
        public int MinFragmentCountForPfGroup { get; set; }
        public abstract IEnumerable<PrecursorFragmentsGroup> PrecursorFragmentGrouping(List<ExtractedIonChromatogram> precursors, List<ExtractedIonChromatogram> fragments);
    }
}
