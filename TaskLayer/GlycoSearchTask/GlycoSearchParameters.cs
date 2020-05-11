using EngineLayer;
using System.Linq;
using UsefulProteomicsDatabases;
using EngineLayer.GlycoSearch;

namespace TaskLayer
{
    public class GlycoSearchParameters
    {
        public GlycoSearchParameters()
        {        
            OGlycanDatabasefile = "OGlycan.gdb";
            NGlycanDatabasefile = "NGlycan182.gdb";
            GlycoSearchType = GlycoSearchType.OGlycanSearch;
            OxoniumIonFilt = true;
            DecoyType = DecoyType.Reverse;
            GlycoSearchTopNum = 50;
            MaximumOGlycanAllowed = 4;
        }

        public string OGlycanDatabasefile { get; set; }
        public string NGlycanDatabasefile { get; set; }
        public GlycoSearchType GlycoSearchType { get; set; }
        public bool OxoniumIonFilt { get; set; }
        public DecoyType DecoyType { get; set; }
        public int GlycoSearchTopNum { get; set; }
        public int MaximumOGlycanAllowed { get; set; }
    }
}