using EngineLayer.CrosslinkSearch;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class XlSearchParameters
    {
        public XlSearchParameters()
        {
            DecoyType = DecoyType.Reverse;
            CrosslinkerType = CrosslinkerType.DSSO;
            RestrictToTopNHits = true;
            CrosslinkSearchTopNum = 300;
            CrosslinkerName = null;
            IsCleavable = false;
            CrosslinkerShortMass = null;
            CrosslinkerLongMass = null;
            CrosslinkerTotalMass = null;
            CrosslinkerDeadEndMassH2O = null;
            CrosslinkerDeadEndMassNH2 = null;
            CrosslinkerDeadEndMassTris = null;
            CrosslinkerResidues = "K";
            CrosslinkerResidues2 = "K";
            XlQuench_H2O = true;
            XlQuench_NH2 = false;
            XlQuench_Tris = true;
            WriteOutputForPercolator = false;
            WritePepXml = true;
            //XlCharge_2_3 = true;
        }
        
        public DecoyType DecoyType { get; set; }
        public CrosslinkerType CrosslinkerType { get; set; }
        public int CrosslinkSearchTopNum { get; set; }
        public string CrosslinkerName { get; set; }
        public double? CrosslinkerTotalMass { get; set; }
        public double? CrosslinkerShortMass { get; set; }
        public double? CrosslinkerLongMass { get; set; }
        public double? CrosslinkerLoopMass { get; set; }
        public string CrosslinkerResidues { get; set; }
        public string CrosslinkerResidues2 { get; set; }
        public double? CrosslinkerDeadEndMassH2O { get; set; }
        public double? CrosslinkerDeadEndMassNH2 { get; set; }
        public double? CrosslinkerDeadEndMassTris { get; set; }
        
        // TODO: 2+ 3+ prime fragments?
        public bool IsCleavable { get; set; }
        public bool RestrictToTopNHits { get; set; }
        public bool WriteOutputForPercolator { get; set; }
        public bool WritePepXml { get; set; }
        public bool XlQuench_H2O { get; set; }
        public bool XlQuench_Tris { get; set; }
        public bool XlQuench_NH2 { get; set; }
        //public bool XlCharge_2_3 { get; set; }
    }
}