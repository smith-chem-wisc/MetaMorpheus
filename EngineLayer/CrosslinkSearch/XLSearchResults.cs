using EngineLayer.CrosslinkSearch;
using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer
{
    public class XLSearchResults : MetaMorpheusEngineResults
    {
        #region Public Constructors

        public XLSearchResults(List<Tuple<PsmCross, PsmCross>> xlpsms, MetaMorpheusEngine searchParams) : base(searchParams)
        {
            XLPsms = xlpsms;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<Tuple<PsmCross, PsmCross>> XLPsms { get; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append(base.ToString());
            return sb.ToString();
        }

        #endregion Public Methods
    }
}