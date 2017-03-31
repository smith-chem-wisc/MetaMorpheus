using System.Collections.Generic;
using System.Text;

namespace EngineLayer.ModernSearch
{
    public class ModernSearchResults : MetaMorpheusEngineResults
    {

        #region Public Constructors

        public ModernSearchResults(List<PsmModern>[] newPsms, ModernSearchEngine s) : base(s)
        {
            this.NewPsms = newPsms;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<PsmModern>[] NewPsms { get; private set; }

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