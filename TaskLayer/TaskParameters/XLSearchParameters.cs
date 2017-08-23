using MzLibUtil;

namespace EngineLayer
{
    public class XlSearchParameters
    {
        public XlSearchParameters()
        {
            SearchDecoy = true;
            CrosslinkerType = CrosslinkerType.DSS;
            CrosslinkSearchTopNum = 50;
            CrosslinkSearchWithAllBeta = false;
            UdXLkerName = null;
            UdXLkerCleavable = false;
            UdXLkerShortMass = null;
            UdXLkerLongMass = null;
            UdXLkerTotalMass = null;
            UdXLkerResidue = 'K';
            XlPrecusorMsTl = new PpmTolerance(10);
            XlBetaPrecusorMsTl = new PpmTolerance(10);
        }

        public bool SearchDecoy { get; set; }
        public bool DoLocalizationAnalysis { get; set; }

        public CrosslinkerType CrosslinkerType { get; set; }
        public int CrosslinkSearchTopNum { get; set; }
        public bool CrosslinkSearchWithAllBeta { get; set; }
        public string UdXLkerName { get; set; }
        public bool UdXLkerCleavable { get; set; }
        public double? UdXLkerTotalMass { get; set; }
        public double? UdXLkerShortMass { get; set; }
        public double? UdXLkerLongMass { get; set; }
        public char UdXLkerResidue { get; set; }
        public Tolerance XlPrecusorMsTl { get; set; }
        public Tolerance XlBetaPrecusorMsTl { get; set; }
        public bool DisposeOfFileWhenDone { get; set; }
    }

}
