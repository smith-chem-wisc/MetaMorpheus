using EngineLayer;
using System.Linq;
using UsefulProteomicsDatabases;
using EngineLayer.CrosslinkSearch;

namespace TaskLayer
{
    public class GlycoSearchParameters
    {
        public GlycoSearchParameters()
        {
            IsOGlycoSearch = true;
            DecoyType = DecoyType.Reverse;
            RestrictToTopNHits = true;
            GlycoSearchTopNum = 50;
            MassDiffAcceptorType = MassDiffAcceptorType.OneMM;
            MaximumOGlycanAllowed = 3;

            AnalyzeOxiniumIons = false;
            FilterScanOxiniumIons = false;
            SearchGlycan182 = false;
        }

        public bool IsOGlycoSearch { get; set; }
        public DecoyType DecoyType { get; set; }
        public int GlycoSearchTopNum { get; set; }
        public bool RestrictToTopNHits { get; set; }
        public MassDiffAcceptorType MassDiffAcceptorType { get; set; }

        public int MaximumOGlycanAllowed { get; set; }

        public bool AnalyzeOxiniumIons { get; set; }
        public bool FilterScanOxiniumIons { get; set; }
        public bool SearchGlycan182 { get; set; }
    }
}