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
            MassDiffAcceptorType = MassDiffAcceptorType.OneMM;
            XlQuench_H2O = true;
            XlQuench_NH2 = false;
            XlQuench_Tris = true;
            WriteOutputForPercolator = false;
            WritePepXml = true;
            AnalyzeOxiniumIons = false;
            FilterScanOxiniumIons = false;
            SearchGlycan182 = false;

        }

        public OpenSearchType OpenSearchType { get; set; }
        public DecoyType DecoyType { get; set; }
        public Crosslinker Crosslinker { get; set; }
        public int CrosslinkSearchTopNum { get; set; }
        public bool RestrictToTopNHits { get; set; }
        public MassDiffAcceptorType MassDiffAcceptorType { get; set; }
        public bool WriteOutputForPercolator { get; set; }
        public bool WritePepXml { get; set; }
        public bool XlQuench_H2O { get; set; }
        public bool XlQuench_Tris { get; set; }
        public bool XlQuench_NH2 { get; set; }
        public bool AnalyzeOxiniumIons { get; set; }
        public bool FilterScanOxiniumIons { get; set; }
        public bool SearchGlycan182 { get; set; }
    }
}