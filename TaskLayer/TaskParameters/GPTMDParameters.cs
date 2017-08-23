using MzLibUtil;
using System.Collections.Generic;
using System;

namespace EngineLayer
{
    public class GPTMDParameters
    {
        public Tolerance PrecursorMassTolerance { get; set; }
        public List<Tuple<string, string>> ListOfModsGptmd { get; set; }
    }
}
