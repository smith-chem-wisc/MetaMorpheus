using System.Collections.Generic;

namespace EngineLayer.Neo
{
    public class NeoPsm
    {
        #region Public Constructors

        public NeoPsm(int scan, double expMass)
        {
            this.scanNumber = scan;
            this.expMass = expMass;
            this.fusionType = FusionCandidate.FusionType.TS; //default
            this.candidates = new List<FusionCandidate>();
        }

        public NeoPsm(int scan, double expMass, InitialID nInfo, InitialID cInfo)
        {
            this.scanNumber = scan;
            this.expMass = expMass;
            this.nInfo = nInfo;
            this.cInfo = cInfo;
            this.candidates = new List<FusionCandidate>();
        }

        #endregion Public Constructors

        #region Public Properties

        public int scanNumber { get; set; }
        public double expMass { get; set; }
        public InitialID nInfo { get; set; }
        public InitialID cInfo { get; set; }
        public List<FusionCandidate> candidates { get; set; }
        public FusionCandidate.FusionType fusionType { get; set; }

        #endregion Public Properties
    }
}