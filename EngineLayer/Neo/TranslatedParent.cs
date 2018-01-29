using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.Neo
{
    public class TranslatedParent : Parent
    {
        public int start { get; set; }
        public int peptideLength { get; set; }
        public FusionCandidate.FusionType translated { get; }

        public TranslatedParent(string id, string seq, int start, int length) : base(id, seq)
        {
            this.start = start;
            this.peptideLength = length;
        }
    }
}
