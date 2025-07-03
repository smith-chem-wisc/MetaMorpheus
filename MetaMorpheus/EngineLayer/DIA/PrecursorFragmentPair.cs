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
        public double Correlation { get; set; }
        public int PrecursorRank { get; set; }
        public int FragmentRank { get; set; }

        public PrecursorFragmentPair(ExtractedIonChromatogram precursorXic, ExtractedIonChromatogram fragmentXic)
        {
            PrecursorXic = precursorXic;
            FragmentXic = fragmentXic;
        }
    }
}

