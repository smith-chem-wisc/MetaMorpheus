namespace EngineLayer.CrosslinkSearch
{
    public enum CrosslinkerType
    {
        DSSO,
        DSS,
        DTSSP,
        BuUrBu,
        DisulfideBond,
        UserDefined
    }

    public class CrosslinkerTypeClass
    {
        #region Public Fields

        public char CrosslinkerModSite;

        #endregion Public Fields

        #region Public Properties

        public string CrosslinkerName { get; set; }
        public bool Cleavable { get; set; }
        public double TotalMass { get; set; }
        public double CleaveMassShort { get; set; }
        public double CleaveMassLong { get; set; }
        public double LoopMass { get; set; }
        public double DeadendMassH2O { get; set; }
        public double DeadendMassNH2 { get; set; }
        public double DeadendMassTris { get; set; }

        #endregion Public Properties

        #region Public Methods

        public CrosslinkerTypeClass SelectCrosslinker(CrosslinkerType name)
        {
            if (name == CrosslinkerType.DSSO)
            {
                CrosslinkerName = "DSSO";
                Cleavable = true;
                TotalMass = 158.0038;
                CleaveMassShort = 54.01056;
                CleaveMassLong = 103.9932;
                CrosslinkerModSite = 'K';
                LoopMass = 159.0012;
                DeadendMassH2O = 176.0143;
                DeadendMassNH2 = 175.0303;
                DeadendMassTris = 279.0777;
                /*Residue.TryGetResidue("K", out CrosslinkerModSite)*/
            }
            if (name == CrosslinkerType.BuUrBu)
            {
                CrosslinkerName = "BuUrBu";
                Cleavable = true;
                TotalMass = 197.0926;
                CleaveMassShort = 111.0320;
                CleaveMassLong = 85.05276;
                CrosslinkerModSite = 'K';
            }
            if (name == CrosslinkerType.DisulfideBond)
            {
                CrosslinkerName = "DisulfideBond";
                Cleavable = true;
                TotalMass = -2.01565;
                CleaveMassShort = -33.98772;
                CleaveMassLong = 31.97207;
                CrosslinkerModSite = 'K';
            }
            if (name == CrosslinkerType.DSS)
            {
                CrosslinkerName = "DSS";
                Cleavable = false;
                TotalMass = 138.06808;
                CrosslinkerModSite = 'K';
                LoopMass = 139.06548;
                DeadendMassH2O = 156.0786;
                DeadendMassNH2 = 155.0946;
                DeadendMassTris = 259.142;
            }
            if (name == CrosslinkerType.DTSSP)
            {
                CrosslinkerName = "DTSSP";
                Cleavable = false;
                TotalMass = 138.06808;
                CrosslinkerModSite = 'K';
            }

            return this;
        }

        #endregion Public Methods
    }
}