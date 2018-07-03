using System.Collections.Generic;

namespace EngineLayer.Neo
{
    public class NeoPsm
    {
        public NeoPsm(int scan, double expMass)
        {
            ScanNumber = scan;
            ExpMass = expMass;
            FusionType = FusionType.TS; //default
            Candidates = new List<FusionCandidate>();
        }

        public NeoPsm(int scan, double expMass, InitialID nInfo, InitialID cInfo)
        {
            ScanNumber = scan;
            ExpMass = expMass;
            NInfo = nInfo;
            CInfo = cInfo;
            Candidates = new List<FusionCandidate>();
        }

        public int ScanNumber { get; set; }
        public double ExpMass { get; set; }
        public InitialID NInfo { get; set; }
        public InitialID CInfo { get; set; }
        public List<FusionCandidate> Candidates { get; set; }
        public FusionType FusionType { get; set; }
    }
}