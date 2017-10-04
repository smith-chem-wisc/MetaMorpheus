using MzLibUtil;
using UsefulProteomicsDatabases;

namespace EngineLayer
{
    public class XlSearchParameters
    {
        #region Public Constructors

        public XlSearchParameters()
        {
            DecoyType = DecoyType.Reverse;
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

        #endregion Public Constructors

        #region Public Properties

        public DecoyType DecoyType { get; set; }
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

        #endregion Public Properties
    }
}