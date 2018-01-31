using Proteomics;
using System.Collections.Generic;

namespace EngineLayer.Neo
{
    //used for parents of a given sequence fragment for FusionCandidate objects
    public class ParentInfo
    {
        #region Public Constructors

        public ParentInfo(List<Protein> proteins, terminal parentType, string seqFound)
        {
            this.theoreticalProteins = proteins;
            this.parentType = parentType;
            this.fragFound = seqFound;
        }

        #endregion Public Constructors

        #region Public Enums

        public enum terminal { N, C };

        #endregion Public Enums

        #region Public Properties

        public List<Protein> theoreticalProteins { get; set; }
        public string fragFound { get; set; }
        public terminal parentType { get; set; }

        #endregion Public Properties

        //What terminus is the fragment from
    }
}