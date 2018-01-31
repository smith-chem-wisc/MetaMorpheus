﻿using System.Collections.Generic;

namespace EngineLayer.Neo
{
    public class CisParent : Parent
    {
        #region Public Constructors

        public CisParent(string id, string seq, List<int> nStart, int nLength, List<int> cStart, int cLength) : base(id, seq)
        {
            this.nStart = nStart;
            this.nLength = nLength;
            this.cStart = cStart;
            this.cLength = cLength;

            //determine cis type
            if (this.nStart[0] <= this.cStart[cStart.Count - 1])
                this.cisType = FusionCandidate.FusionType.NC;
            else
                this.cisType = FusionCandidate.FusionType.RC;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<int> nStart { get; set; }
        public int nLength { get; set; }
        public List<int> cStart { get; set; }
        public int cLength { get; set; }
        public FusionCandidate.FusionType cisType { get; }

        #endregion Public Properties
    }
}