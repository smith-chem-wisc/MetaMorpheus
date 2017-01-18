using System.Collections.Generic;
using System.Text;

namespace InternalLogicEngineLayer
{
    public class ModernSearchResults : MyResults
    {
        #region Public Constructors

        public ModernSearchResults(List<ModernSpectrumMatch>[] newPsms, ModernSearchEngine s) : base(s)
        {
            this.NewPsms = newPsms;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<ModernSpectrumMatch>[] NewPsms { get; private set; }

        #endregion Public Properties

        #region Protected Methods

        protected override string StringForOutput
        {
            get
            {
                var sb = new StringBuilder();
                return sb.ToString();
            }
        }

        #endregion Protected Methods
    }
}