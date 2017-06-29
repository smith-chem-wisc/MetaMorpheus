using System.Text;

namespace EngineLayer.ModernSearch
{
    public class ModernSearchResults : MetaMorpheusEngineResults
    {

        #region Public Constructors

        public ModernSearchResults(PsmParent[][] newPsms, ModernSearchEngine s) : base(s)
        {
            this.OuterPsms = newPsms;
        }

        #endregion Public Constructors

        #region Public Properties

        public PsmParent[][] OuterPsms { get; }

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