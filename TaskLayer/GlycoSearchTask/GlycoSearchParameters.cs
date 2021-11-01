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
            MixedGlycoAllowed = false;
            OxoniumIonFilt = true;
            DecoyType = DecoyType.Reverse;
            GlycoSearchTopNum = 20;
            MaximumOGlycanAllowed = 4;
            MaximumNGlycanAllowed = 1;
            Indexing_by_ion = true;
            IndexingChildScan = false;
            IndexingChildScanDiffIndex = false;
        }

        public string OGlycanDatabasefile { get; set; }
        public string NGlycanDatabasefile { get; set; }
        public GlycoSearchType GlycoSearchType { get; set; }

        public bool MixedGlycoAllowed { get; set; } 
        public bool OxoniumIonFilt { get; set; }
        public DecoyType DecoyType { get; set; }
        public int GlycoSearchTopNum { get; set; }
        public int MaximumOGlycanAllowed { get; set; }

        public int MaximumNGlycanAllowed { get; set; }

        public bool Indexing_by_ion { get; set; }
        public bool IndexingChildScan { get; set; }
        public bool IndexingChildScanDiffIndex { get; set; }
    }
}