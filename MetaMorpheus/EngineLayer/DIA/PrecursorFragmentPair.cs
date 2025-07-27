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

    }
}

