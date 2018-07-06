using Proteomics;
using System.Collections.Generic;

namespace EngineLayer.Neo
{
    //used for parents of a given sequence fragment for FusionCandidate objects
    public class ParentInfo
    {
        public ParentInfo(List<Protein> proteins, terminal parentType, string seqFound)
        {
            this.theoreticalProteins = proteins;
            this.parentType = parentType;
            this.fragFound = seqFound;
        }

        public enum terminal { N, C };

        public List<Protein> theoreticalProteins { get; set; }
        public string fragFound { get; set; }
        public terminal parentType { get; set; }

        //What terminus is the fragment from
    }
}