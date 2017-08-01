using System.Collections.Generic;
using System.Text;
using System;

namespace EngineLayer.CrosslinkSearch
{
    public class CrosslinkSearchResults : MetaMorpheusEngineResults
    {
        #region Public Constructors

        public CrosslinkSearchResults(List<Tuple<PsmCross, PsmCross>> newPsms, CrosslinkSearchEngine s) : base(s)
        {
            this.NewPsms = newPsms;
        }

        public CrosslinkSearchResults(List<Tuple<PsmCross, PsmCross>> newPsms, CrosslinkSearchEngine2 s) : base(s)
        {
            this.NewPsms = newPsms;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<Tuple<PsmCross, PsmCross>> NewPsms { get; private set; }

        #endregion Public Properties
    }
}
