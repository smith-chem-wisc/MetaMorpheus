using System.Text;

namespace EngineLayer.ClassicSearch
{
    public class ClassicSearchResults : MetaMorpheusEngineResults
    {

        #region Public Constructors

        public ClassicSearchResults(PsmParent[][] outerPsms, ClassicSearchEngine searchParams) : base(searchParams)
        {
            OuterPsms = outerPsms;
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