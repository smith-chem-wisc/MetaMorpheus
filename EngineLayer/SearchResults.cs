using System.Text;

namespace EngineLayer
{
    public class SearchResults : MetaMorpheusEngineResults
    {
        #region Public Constructors

        public SearchResults(Psm[][] psms, MetaMorpheusEngine searchParams) : base(searchParams)
        {
            Psms = psms;
        }

        #endregion Public Constructors

        #region Public Properties

        public Psm[][] Psms { get; }

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