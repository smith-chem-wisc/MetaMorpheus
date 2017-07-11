namespace EngineLayer
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
        public string CrosslinkerName { get; set; }
        public bool Cleavable { get; set; }
        public double TotalMass { get; set; }
        public double CleaveMassShort { get; set; }
        public double CleaveMassLong { get; set; }
        public string CrosslinkerModSite;

        public CrosslinkerTypeClass SelectCrosslinker(CrosslinkerType name)
        {

            if (name == CrosslinkerType.DSSO)
            {
                CrosslinkerName = "DSSO";
                Cleavable = true;
                TotalMass = 158.0038;
                CleaveMassShort = 54.01056;
                CleaveMassLong = 103.9932;
                CrosslinkerModSite = "K";
                /*Residue.TryGetResidue("K", out CrosslinkerModSite)*/
                ;
            }
            //if (name == CrosslinkerType.BuUrBu)
            //{
            //    CrosslinkerName = "BuUrBu";
            //    Cleavable = true;
            //    TotalMass = 197.0926;
            //    CleaveMassShort = 111.0320;
            //    CleaveMassLong = 85.05276;
            //    Residue.TryGetResidue("K", out CrosslinkerModSite);
            //}
            //if (name == CrosslinkerType.DisulfideBond)
            //{
            //    CrosslinkerName = "DisulfideBond";
            //    Cleavable = true;
            //    TotalMass = -2.01565;
            //    CleaveMassShort = -33.98772;
            //    CleaveMassLong = 31.97207;
            //    Residue.TryGetResidue("K", out CrosslinkerModSite);
            //}
            if (name == CrosslinkerType.DSS)
            {
                CrosslinkerName = "DSS";
                Cleavable = false;
                TotalMass = 138.06808;
                CrosslinkerModSite = "K";
                //Residue.TryGetResidue("K", out CrosslinkerModSite);
            }
            //if (name == CrosslinkerType.DTSSP)
            //{
            //    CrosslinkerName = "DTSSP";
            //    Cleavable = false;
            //    TotalMass = 138.06808;
            //    Residue.TryGetResidue("K", out CrosslinkerModSite);
            //}

            return this;
        }
    }
}
