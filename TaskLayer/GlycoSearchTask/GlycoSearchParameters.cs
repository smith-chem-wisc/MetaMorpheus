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
            GlycanDatabasefile = "OGlycan.gdb";
            IsOGlycoSearch = true;
            OxiniumIonFilt = true;
            DecoyType = DecoyType.Reverse;
            GlycoSearchTopNum = 50;
            MaximumOGlycanAllowed = 4;
        }

        public string GlycanDatabasefile { get; set; }
        public bool IsOGlycoSearch { get; set; }
        public bool OxiniumIonFilt { get; set; }
        public DecoyType DecoyType { get; set; }
        public int GlycoSearchTopNum { get; set; }
        public int MaximumOGlycanAllowed { get; set; }
    }
}