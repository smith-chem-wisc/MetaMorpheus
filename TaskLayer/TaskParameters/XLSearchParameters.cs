using MzLibUtil;

namespace EngineLayer
{
    public class XlSearchParameters
    {
        #region Public Constructors

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

            XlOutAll = true;
            XlOutCrosslink = true;
            XlOutCLMSVault = false;
            XlOutPercolator = false;
            XlOutPepXML = false;
        }

        #endregion Public Constructors

        #region Public Properties

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

        public bool XlOutPercolator { get; set; }
        public bool XlOutCLMSVault { get; set; }
        public bool XlOutAll { get; set; }
        public bool XlOutPepXML { get; set; }
        public bool XlOutCrosslink { get; set; }
        #endregion Public Properties
    }
}