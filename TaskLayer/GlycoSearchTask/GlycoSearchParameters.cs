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
            GlycanDatabasefile = "OGlycan.gdb";

            DecoyType = DecoyType.Reverse;
            GlycoSearchTopNum = 50;

            MaximumOGlycanAllowed = 3;

            AnalyzeOxiniumIons = false;
            FilterScanOxiniumIons = false;
        }

        public string GlycanDatabasefile { get; set; }
        public bool IsOGlycoSearch { get; set; }
        public DecoyType DecoyType { get; set; }
        public int GlycoSearchTopNum { get; set; }


        public int MaximumOGlycanAllowed { get; set; }

        public bool AnalyzeOxiniumIons { get; set; }
        public bool FilterScanOxiniumIons { get; set; }
    }
}