using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using MzLibUtil;


namespace EngineLayer.DIA
{
    public class PrecursorFragmentPair
    {
        public ExtractedIonChromatogram PrecursorXic { get; set; }
        public ExtractedIonChromatogram FragmentXic { get; set; }
        public double? Correlation { get; set; }
        public double? Overlap { get; set; } 
        public int? PrecursorRank { get; set; }
        public int? FragmentRank { get; set; }

        public PrecursorFragmentPair(ExtractedIonChromatogram precursorXic, ExtractedIonChromatogram fragmentXic, double? correlation = null, double? overlap = null)
        {
            PrecursorXic = precursorXic;
            FragmentXic = fragmentXic;
            Correlation = correlation;
            Overlap = overlap;
        }

        public static void SetPrecursorRankForPfPairs(List<PrecursorFragmentPair> allPfPairs)
        {
            var fragmentPairMap = new Dictionary<ExtractedIonChromatogram, List<PrecursorFragmentPair>>();
            foreach (var pair in allPfPairs)
            {
                if (!fragmentPairMap.ContainsKey(pair.FragmentXic))
                {
                    fragmentPairMap[pair.FragmentXic] = new List<PrecursorFragmentPair>();
                }
                fragmentPairMap[pair.FragmentXic].Add(pair);
            }
            foreach (var fragmentXic in fragmentPairMap.Keys)
            {
                fragmentPairMap[fragmentXic].MaxBy(pf => pf.Correlation);
                for (int i = 0; i < fragmentPairMap[fragmentXic].Count; i++)
                {
                    fragmentPairMap[fragmentXic][i].PrecursorRank = i + 1;
                }
            }
        }
    }
}

