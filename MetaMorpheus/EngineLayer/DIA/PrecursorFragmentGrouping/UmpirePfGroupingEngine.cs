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
    /// <summary>
    /// UmpirePfGroupingEngine uses algorithms from DIA-Umpire for grouping XICs where overlap ratio and Pearson's correlation are calculated differently than in XicGroupingEngine. 
    /// Reference to DIA-Umpire: Tsou et al.,Nat Methods. 2015 Mar;12(3):258-64. doi: 10.1038/nmeth.3255. Epub 2015 Feb 2. PMID: 25665051; PMCID: PMC4346634.
    /// </summary>
    public class UmpirePfGroupingEngine : XicGroupingEngine
    {
        public int NoPointPerInterval { get; set; }
        public UmpirePfGroupingEngine(int noPointPerInterval, float apexRTTolerance, double overlapThreshold, double correlationThreshold, int maxThreadsForGrouping = 1, int minFragmentCountForGrouping = 0, int? precursorRankThreshold = null, int? fragmentRankThreshold = null)
            : base(apexRTTolerance, overlapThreshold, correlationThreshold, maxThreadsForGrouping, minFragmentCountForGrouping, precursorRankThreshold, fragmentRankThreshold)
        {
            NoPointPerInterval = noPointPerInterval;
        }

        /// <summary>
        /// The general workflow of grouping is the same as XicGroupingEngine but the overlap ratio and correlation calculations use DIA-Umpire's algorithms.
        /// </summary>
        public override PrecursorFragmentsGroup GroupFragmentsForOnePrecursor(ExtractedIonChromatogram precursorXic, IEnumerable<ExtractedIonChromatogram> fragmentXics)
        {
            var pfPairs = new List<PrecursorFragmentPair>();
            foreach (var fragmentXic in fragmentXics)
            {
                //umpire style overlap and correlation calculation
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
