using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.CrosslinkSearch
{
    public class CrosslinkInfo
    {
        public double XLTotalScore { get; set; } //alpha + beta psmCross
        public double XLQvalueTotalScore { get; set; } //Calc based on XLtotalScore for Qvalue
        public int XlProteinPos { get; set; }
        public int[] XlRank { get; set; } //only contain 2 intger, consider change to Tuple
        public string ParentIonExist { get; set; }
        public int ParentIonExistNum { get; set; }
        public List<int> ParentIonMaxIntensityRanks { get; set; }
        public PsmCrossType CrossType { get; set; }
    }
}
