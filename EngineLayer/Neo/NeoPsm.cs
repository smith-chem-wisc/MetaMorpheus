using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.Neo
{
    public class NeoPsm
    {
        public string file { get; set; }
        public int scanNumber { get; set; }
        public double expMass { get; set; }
        public InitialID nInfo { get; set; }
        public InitialID cInfo { get; set; }
        public List<FusionCandidate> candidates { get; set; }
        public FusionCandidate.FusionType fusionType { get; set; }


        public NeoPsm(string file, int scan, double expMass)
        {
            this.file = file;
            this.scanNumber = scan;
            this.expMass = expMass;
            this.fusionType = FusionCandidate.FusionType.TS; //default
            this.candidates = new List<FusionCandidate>();
        }

        public NeoPsm(string file, int scan, double expMass, InitialID nInfo, InitialID cInfo)
        {
            this.file = file;
            this.scanNumber = scan;
            this.expMass = expMass;
            this.nInfo = nInfo;
            this.cInfo = cInfo;
            this.candidates = new List<FusionCandidate>();
        }
    }
}
