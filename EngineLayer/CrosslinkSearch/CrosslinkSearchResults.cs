using System;
using System.Collections.Generic;

namespace EngineLayer.CrosslinkSearch
{
    public class CrosslinkSearchResults : MetaMorpheusEngineResults
    {
        #region Public Constructors

        public CrosslinkSearchResults(List<PsmCross> newPsms, CrosslinkSearchEngine s) : base(s)
        {
            this.NewPsms = newPsms;
        }

        public CrosslinkSearchResults(List<PsmCross> newPsms, CrosslinkSearchEngine2 s) : base(s)
        {
            this.NewPsms = newPsms;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<PsmCross> NewPsms { get; private set; }

        #endregion Public Properties
    }
}