using MzLibUtil;
using System.Collections.Generic;
using System;

namespace EngineLayer
{
    //XLSearch Parameters
    public class XLSearchParameters
    {
        public XLSearchParameters()
        {
            DoParsimony = false;
            NoOneHitWonders = false;
            ModPeptidesAreUnique = false;
            DoQuantification = false;
            QuantifyPpmTol = 5;
            SearchDecoy = true;
            DoHistogramAnalysis = false;

            DoLocalizationAnalysis = true;
            WritePrunedDatabase = false;
            KeepAllUniprotMods = true;
            CrosslinkerType = CrosslinkerType.DSS;
            CrosslinkSearchTopNum = 50;
            CrosslinkSearchWithAllBeta = false;
            UdXLkerName = null;
            UdXLkerCleavable = false;
            UdXLkerShortMass = null;
            UdXLkerLongMass = null;
            UdXLkerTotalMass = null;
            UdXLkerResidue = 'K';
            XLprecusorMsTl = new PpmTolerance(10);
            XLBetaPrecusorMsTl = new PpmTolerance(10);
            DisposeOfFileWhenDone = true;
        }
        public bool DoParsimony { get; set; }
        public bool ModPeptidesAreUnique { get; set; }
        public bool NoOneHitWonders { get; set; }
        public bool MatchBetweenRuns { get; set; }
        public double QuantifyPpmTol { get; set; }
        public bool DoHistogramAnalysis { get; set; }
        public bool SearchDecoy { get; set; }

        public bool WritePrunedDatabase { get; set; }
        public bool KeepAllUniprotMods { get; set; }

        public bool DoLocalizationAnalysis { get; set; }
        public bool DoQuantification { get; set; }

        public CrosslinkerType CrosslinkerType { get; set; }
        public int CrosslinkSearchTopNum { get; set; }
        public bool CrosslinkSearchWithAllBeta { get; set; }
        public string UdXLkerName { get; set; }
        public bool UdXLkerCleavable { get; set; }
        public double? UdXLkerTotalMass { get; set; }
        public double? UdXLkerShortMass { get; set; }
        public double? UdXLkerLongMass { get; set; }
        public char UdXLkerResidue { get; set; }
        public Tolerance XLprecusorMsTl { get; set; }
        public Tolerance XLBetaPrecusorMsTl { get; set; }
        public bool DisposeOfFileWhenDone { get; set; }
    }

}
