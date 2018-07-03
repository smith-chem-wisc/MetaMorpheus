using Proteomics;
using System.Collections.Generic;

namespace EngineLayer.Neo
{
    /// <summary>
    /// Used for parents of a given sequence fragment for FusionCandidate objects
    /// </summary>
    public class ParentInfo
    {
        public ParentInfo(List<Protein> proteins, Terminal parentType, string seqFound)
        {
            TheoreticalProteins = proteins;
            ParentType = parentType;
            FragFound = seqFound;
        }

        public enum Terminal { N, C };

        public List<Protein> TheoreticalProteins { get; set; }
        public string FragFound { get; set; }
        public Terminal ParentType { get; set; }

        //What terminus is the fragment from
    }
}