using MassSpectrometry;


namespace EngineLayer.DIA
{
    /// <summary>
    /// PrecursorFragmentPair class represents a pair of precursor XIC and fragment XIC, along with their relationship metrics such as correlation and overlap.
    /// <summary>
    public class PrecursorFragmentPair
    {
        public ExtractedIonChromatogram PrecursorXic { get; set; }
        public ExtractedIonChromatogram FragmentXic { get; set; }
        public double? Correlation { get; set; }
        public double? Overlap { get; set; } 
        public PrecursorFragmentPair(ExtractedIonChromatogram precursorXic, ExtractedIonChromatogram fragmentXic, double? correlation = null, double? overlap = null)
        {
            PrecursorXic = precursorXic;
            FragmentXic = fragmentXic;
            Correlation = correlation;
            Overlap = overlap;
        }
    }
}

