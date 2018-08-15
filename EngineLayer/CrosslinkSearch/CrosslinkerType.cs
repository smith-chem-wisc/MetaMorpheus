namespace EngineLayer.CrosslinkSearch
{
    public enum CrosslinkerType
    {
        DSSO,
        DSS,
        DisulfideBond,
        UserDefined
    }

    public class CrosslinkerTypeClass
    {
        public CrosslinkerTypeClass(string crosslinkerModSites, string crosslinkerModSites2, string crosslinkerName, bool cleavable, double totalMass,
            double cleaveMassShort, double cleaveMassLong, double loopMass, double deadendMassH2O, double deadendMassNH2, double deadendMassTris)
        {
            CrosslinkerModSites = crosslinkerModSites;
            CrosslinkerModSites2 = crosslinkerModSites2;
            CrosslinkerName = crosslinkerName;
            Cleavable = cleavable;
            TotalMass = totalMass;
            CleaveMassShort = cleaveMassShort;
            CleaveMassLong = cleaveMassLong;
            LoopMass = loopMass;
            DeadendMassH2O = deadendMassH2O;
            DeadendMassNH2 = deadendMassNH2;
            DeadendMassTris = deadendMassTris;
        }

        public CrosslinkerTypeClass()
        {
        }

        public string CrosslinkerModSites { get; set; }
        public string CrosslinkerModSites2 { get; set; }
        public string CrosslinkerName { get; set; }
        public bool Cleavable { get; set; }
        public double TotalMass { get; set; }
        public double CleaveMassShort { get; set; }
        public double CleaveMassLong { get; set; }
        public double LoopMass { get; set; }
        public double DeadendMassH2O { get; set; }
        public double DeadendMassNH2 { get; set; }
        public double DeadendMassTris { get; set; }

        public CrosslinkerTypeClass SelectCrosslinker(CrosslinkerType name)
        {
            if (name == CrosslinkerType.DSSO)
            {
                CrosslinkerName = "DSSO";
                Cleavable = true;
                TotalMass = 158.0038;
                CleaveMassShort = 54.01056;
                CleaveMassLong = 103.9932;
                CrosslinkerModSites = "K";
                CrosslinkerModSites2 = "K";
                LoopMass = 158.0038;
                DeadendMassH2O = 176.0143;
                DeadendMassNH2 = 175.0303;
                DeadendMassTris = 279.0777;
                /*Residue.TryGetResidue("K", out CrosslinkerModSite)*/
            }
            if (name == CrosslinkerType.DisulfideBond)
            {
                CrosslinkerName = "DisulfideBond";
                Cleavable = true;
                TotalMass = -2.01565;
                CleaveMassShort = -33.98772;
                CleaveMassLong = 31.97207;
                CrosslinkerModSites = "C";
                CrosslinkerModSites2 = "C";
            }
            if (name == CrosslinkerType.DSS)
            {
                CrosslinkerName = "DSS";
                Cleavable = false;
                TotalMass = 138.06808;
                CrosslinkerModSites = "K";
                CrosslinkerModSites2 = "K";
                LoopMass = 138.06808;
                DeadendMassH2O = 156.0786;
                DeadendMassNH2 = 155.0946;
                DeadendMassTris = 259.142;
            }

            return this;
        }
    }
}