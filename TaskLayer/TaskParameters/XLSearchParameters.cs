using EngineLayer;
using MzLibUtil;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class XlSearchParameters
    {
        #region Public Constructors

        public XlSearchParameters()
        {
            DecoyType = DecoyType.Reverse;
            CrosslinkerType = CrosslinkerType.DSSO;
            CrosslinkSearchTop = false;
            CrosslinkSearchTopNum = 300;
            CrosslinkSearchWithAllBeta = false;
            UdXLkerName = null;
            UdXLkerCleavable = false;
            UdXLkerShortMass = null;
            UdXLkerLongMass = null;
            UdXLkerTotalMass = null;
            UdXLkerResidue = 'K';
            XlPrecusorMsTl = new PpmTolerance(10);
            XlBetaPrecusorMsTl = new PpmTolerance(10);
            XlQuench_H2O = true;
            XlQuench_NH2 = false;
            XlQuench_Tris = true;

            XlOutAll = true;
            XlOutCrosslink = true;
            XlOutPercolator = false;
            XlOutPepXML = false;

            XlCharge_2_3 = true;
            XlCharge_2_3_PrimeFragment = false;
        }

        #endregion Public Constructors

        #region Public Properties

        public DecoyType DecoyType { get; set; }
        public bool DoLocalizationAnalysis { get; set; }

        public CrosslinkerType CrosslinkerType { get; set; }
        public bool CrosslinkSearchTop { get; set; }
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
        public bool XlOutAll { get; set; }
        public bool XlOutPepXML { get; set; }
        public bool XlOutCrosslink { get; set; }

        public bool XlQuench_H2O { get; set; }
        public bool XlQuench_Tris { get; set; }
        public bool XlQuench_NH2 { get; set; }

        public bool XlCharge_2_3 { get; set; }
        public bool XlCharge_2_3_PrimeFragment { get; set; }

        #endregion Public Properties
    }
}