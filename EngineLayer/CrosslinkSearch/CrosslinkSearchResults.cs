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

        //#region Public Methods

        //public override string ToString()
        //{
        //    var sb = new StringBuilder();
        //    sb.Append(base.ToString());
        //    return sb.ToString();
        //}

        //#endregion Public Methods
    }
}
