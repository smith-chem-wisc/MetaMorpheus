using EngineLayer.CrosslinkSearch;
using MzLibUtil;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class XlSearchParameters
    {
        #region Public Constructors

        public XlSearchParameters()
        {
            FragmentationType = FragmentaionType.MS2_HCD;
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
            UdXLkerLoopMass = null;
            UdXLkerResidues = "K";
            UdXLkerResidues2 = "K";
            UdXLkerDeadendMassH2O = null;
            UdXLkerDeadendMassNH2 = null;
            UdXLkerDeadendMassTris = null;
            XlPrecusorMsTl = new PpmTolerance(10);
            XlBetaPrecusorMsTl = new PpmTolerance(10);
            XlQuench_H2O = true;
            XlQuench_NH2 = false;
            XlQuench_Tris = true;

            XlOutAll = false;
            XlOutCrosslink = true;
            XlOutPercolator = false;
            XlOutPepXML = false;

            XlCharge_2_3 = true;
            XlCharge_2_3_PrimeFragment = false;
        }

        #endregion Public Constructors

        #region Public Properties

        public FragmentaionType FragmentationType { get; set; }
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
        public double? UdXLkerLoopMass { get; set; }
        public string UdXLkerResidues { get; set; }
        public string UdXLkerResidues2 { get; set; }
        public double? UdXLkerDeadendMassH2O { get; set; }
        public double? UdXLkerDeadendMassNH2 { get; set; }
        public double? UdXLkerDeadendMassTris { get; set; }
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