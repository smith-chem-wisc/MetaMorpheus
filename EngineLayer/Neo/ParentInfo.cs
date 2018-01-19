using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.Neo
{
    //used for parents of a given sequence fragment for FusionCandidate objects
    public class ParentInfo
    {
        public List<Parent> theoreticalProteins { get; set; }
        public string fragFound { get; set; }
        public enum terminal { N, C };
        public terminal parentType { get; set; } //What terminus is the fragment from
        public ParentInfo(List<Parent> proteins, terminal parentType, string seqFound)
        {
            this.theoreticalProteins = proteins;
            this.parentType = parentType;
            this.fragFound = seqFound;
        }
    }
}
