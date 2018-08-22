using EngineLayer.CrosslinkSearch;
using MzLibUtil;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class XlSearchParameters
    {
        public XlSearchParameters()
        {
            SearchGlyco = false;
            SearchGlycoWithBgYgIndex = false;
            DecoyType = DecoyType.Reverse;
            CrosslinkerType = CrosslinkerType.DSSO;
            CrosslinkSearchTop = true;
            CrosslinkSearchTopNum = 300;
            UdXLkerName = null;
            UdXLkerCleavable = false;
            UdXLkerShortMass = null;
            UdXLkerLongMass = null;
            UdXLkerTotalMass = null;
            UdXLkerResidues = "K";
            UdXLkerResidues2 = "K";
            XlPrecusorMsTl = new PpmTolerance(10);
            XlQuench_H2O = true;
            XlQuench_NH2 = false;
            XlQuench_Tris = true;

            XlOutAll = false;
            XlOutCrosslink = true;
            XlOutPercolator = false;
            XlOutPepXML = false;

            XlCharge_2_3 = true;
        }

        public bool SearchGlyco { get; set; }
        public bool SearchGlycoWithBgYgIndex { get ;set;}
        public DecoyType DecoyType { get; set; }
        public CrosslinkerType CrosslinkerType { get; set; }
        public bool CrosslinkSearchTop { get; set; }
        public int CrosslinkSearchTopNum { get; set; }
        public string UdXLkerName { get; set; }
        public bool UdXLkerCleavable { get; set; }
        public double? UdXLkerTotalMass { get; set; }
        public double? UdXLkerShortMass { get; set; }
        public double? UdXLkerLongMass { get; set; }
        public string UdXLkerResidues { get; set; }
        public string UdXLkerResidues2 { get; set; }

        public Tolerance XlPrecusorMsTl { get; set; }
        public bool DisposeOfFileWhenDone { get; set; }

        public bool XlOutPercolator { get; set; }
        public bool XlOutAll { get; set; }
        public bool XlOutPepXML { get; set; }
        public bool XlOutCrosslink { get; set; }

        public bool XlQuench_H2O { get; set; }
        public bool XlQuench_Tris { get; set; }
        public bool XlQuench_NH2 { get; set; }

        public bool XlCharge_2_3 { get; set; }
    }
}