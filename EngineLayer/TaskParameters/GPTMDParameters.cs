using MzLibUtil;
using System.Collections.Generic;
using System;

namespace EngineLayer
{
    //GPTMD Parameters
    public class GPTMDParameters
    {
        public GPTMDParameters()
        {
            PrecursorMassTolerance = new PpmTolerance(2);
            //    ListOfModsGptmd = GlobalTaskLevelSettings.AllModsKnown.Where(b =>
            //        b.modificationType.Equals("Glycan") ||
            //        b.modificationType.Equals("Mod") ||
            //        b.modificationType.Equals("PeptideTermMod") ||
            //        b.modificationType.Equals("Metal") ||
            //        b.modificationType.Equals("ProteinTermMod")).Select(b => new Tuple<string, string>(b.modificationType, b.id)).ToList();
        }
        public Tolerance PrecursorMassTolerance { get; set; }
        public List<Tuple<string, string>> ListOfModsGptmd { get; set; }
    }
}
