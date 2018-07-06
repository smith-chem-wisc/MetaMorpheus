using EngineLayer.CrosslinkSearch;
using MzLibUtil;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class XlSearchParameters
    {
        public FragmentaionType FragmentationType { get; set; } = FragmentaionType.MS2_HCD;
        public DecoyType DecoyType { get; set; } = DecoyType.Reverse;
        public bool DoLocalizationAnalysis { get; set; }

        public CrosslinkerType CrosslinkerType { get; set; } = CrosslinkerType.DSSO;
        public bool CrosslinkSearchTop { get; set; }
        public int CrosslinkSearchTopNum { get; set; } = 300;
        public bool CrosslinkSearchWithAllBeta { get; set; }
        public string UdXLkerName { get; set; }
        public bool UdXLkerCleavable { get; set; }
        public double? UdXLkerTotalMass { get; set; }
        public double? UdXLkerShortMass { get; set; }
        public double? UdXLkerLongMass { get; set; }
        public double? UdXLkerLoopMass { get; set; }
        public string UdXLkerResidues { get; set; } = "K";
        public string UdXLkerResidues2 { get; set; } = "K";
        public double? UdXLkerDeadendMassH2O { get; set; }
        public double? UdXLkerDeadendMassNH2 { get; set; }
        public double? UdXLkerDeadendMassTris { get; set; }
        public Tolerance XlPrecusorMsTl { get; set; } = new PpmTolerance(10);
        public Tolerance XlBetaPrecusorMsTl { get; set; } = new PpmTolerance(10);
        public bool DisposeOfFileWhenDone { get; set; }

        public bool XlOutPercolator { get; set; }
        public bool XlOutAll { get; set; }
        public bool XlOutPepXML { get; set; }
        public bool XlOutCrosslink { get; set; } = true;

        public bool XlQuench_H2O { get; set; } = true;
        public bool XlQuench_Tris { get; set; } = true;
        public bool XlQuench_NH2 { get; set; }

        public bool XlCharge_2_3 { get; set; } = true;
        public bool XlCharge_2_3_PrimeFragment { get; set; }
    }
}