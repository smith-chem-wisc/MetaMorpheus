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
            CrosslinkSearchTopNum = 50;
            CrosslinkAtCleavageSite = false;
            XlQuench_H2O = true;
            XlQuench_NH2 = false;
            XlQuench_Tris = false;
            WritePepXml = true;
            WriteOutputForPercolator = false;
        }
        
        public DecoyType DecoyType { get; set; }
        public Crosslinker Crosslinker { get; set; }
        public int CrosslinkSearchTopNum { get; set; }
        public bool CrosslinkAtCleavageSite { get; set; }
        public bool WritePepXml { get; set; }
        public bool WriteOutputForPercolator { get; set; }
        public bool XlQuench_H2O { get; set; }
        public bool XlQuench_Tris { get; set; }
        public bool XlQuench_NH2 { get; set; }
    }
}