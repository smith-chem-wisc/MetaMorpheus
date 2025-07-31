using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public abstract class PfGroupingEngine
    {
        /// <summary>
        /// PfGroupingEngine defines the process of grouping precursor and fragment XICs into PrecursorFragmentsGroup objects. 
        /// It should have a method that returns all PrecursorFragmentsGroup objects that can be found in a given set of precursor and fragment XICs.
        /// <summary>
        public int MaxThreadsForGrouping { get; set; }
        public int MinFragmentCountForPfGroup { get; set; }
        public abstract List<PrecursorFragmentsGroup> PrecursorFragmentGrouping(List<ExtractedIonChromatogram> precursors, List<ExtractedIonChromatogram> fragments);
    }
}
