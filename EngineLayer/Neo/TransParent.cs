using System.Collections.Generic;

namespace EngineLayer.Neo
{
    public class TransParent : Parent
    {
        public TransParent(string id, string seq, List<int> start, int length, ParentInfo.Terminal terminal)
            : base(id, seq)
        {
            Start = start;
            PeptideLength = length;
            Terminal = terminal;
        }

        public List<int> Start { get; set; }
        public int PeptideLength { get; set; }
        public ParentInfo.Terminal Terminal { get; set; }
    }
}