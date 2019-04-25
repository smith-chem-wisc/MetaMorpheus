using EngineLayer;
using System.Linq;
using UsefulProteomicsDatabases;
using EngineLayer.CrosslinkSearch;

namespace TaskLayer
{
    public class XlSearchParameters
    {
        public XlSearchParameters()
        {
            OpenSearchType = OpenSearchType.Crosslink;
            DecoyType = DecoyType.Reverse;
            Crosslinker = GlobalVariables.Crosslinkers.First();
            RestrictToTopNHits = true;
            CrosslinkSearchTopNum = 300;
            XlQuench_H2O = true;
            XlQuench_NH2 = false;
            XlQuench_Tris = true;
            WriteOutputForPercolator = false;
            WritePepXml = true;
            OnlyAnalyzeOxiniumIons = false;
        }

        public OpenSearchType OpenSearchType { get; set; }
        public DecoyType DecoyType { get; set; }
        public Crosslinker Crosslinker { get; set; }
        public int CrosslinkSearchTopNum { get; set; }
        public bool RestrictToTopNHits { get; set; }
        public bool WriteOutputForPercolator { get; set; }
        public bool WritePepXml { get; set; }
        public bool XlQuench_H2O { get; set; }
        public bool XlQuench_Tris { get; set; }
        public bool XlQuench_NH2 { get; set; }
        public bool OnlyAnalyzeOxiniumIons { get; set; }
    }
}