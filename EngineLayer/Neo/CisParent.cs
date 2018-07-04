using System.Collections.Generic;

namespace EngineLayer.Neo
{
    public class CisParent : Parent
    {
        public CisParent(string id, string seq, List<int> nStart, int nLength, List<int> cStart, int cLength) : base(id, seq)
        {
            NStart = nStart;
            NLength = nLength;
            CStart = cStart;
            CLength = cLength;

            //determine cis type
            CisType = NStart[0] <= CStart[cStart.Count - 1] ? FusionType.NC : FusionType.RC;
        }

        public List<int> NStart { get; set; }
        public int NLength { get; set; }
        public List<int> CStart { get; set; }
        public int CLength { get; set; }
        public FusionType CisType { get; }
    }
}