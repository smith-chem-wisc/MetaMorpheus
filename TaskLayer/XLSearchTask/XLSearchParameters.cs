using EngineLayer;
using System.Linq;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class XlSearchParameters
    {
        public XlSearchParameters()
        {
            DecoyType = DecoyType.Reverse;
            Crosslinker = GlobalVariables.Crosslinkers.First();
            RestrictToTopNHits = true;
            CrosslinkSearchTopNum = 300;
            CrosslinkAtCleavageSite = false;
            XlQuench_H2O = true;
            XlQuench_NH2 = false;
            XlQuench_Tris = true;
            WriteOutputForPercolator = false;
            WritePepXml = true;
        }
        
        public DecoyType DecoyType { get; set; }
        public Crosslinker Crosslinker { get; set; }
        public int CrosslinkSearchTopNum { get; set; }
        public bool RestrictToTopNHits { get; set; }
        public bool CrosslinkAtCleavageSite { get; set; }
        public bool WriteOutputForPercolator { get; set; }
        public bool WritePepXml { get; set; }
        public bool XlQuench_H2O { get; set; }
        public bool XlQuench_Tris { get; set; }
        public bool XlQuench_NH2 { get; set; }
    }
}