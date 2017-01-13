using System.Text;

namespace InternalLogicEngineLayer
{
    public class ClassicSearchResults : MyResults
    {
        #region Internal Constructors

        internal ClassicSearchResults(ClassicSearchEngine searchParams) : base(searchParams)
        {
        }

        #endregion Internal Constructors

        #region Public Properties

        public ClassicSpectrumMatch[][] outerPsms { get; set; }

        #endregion Public Properties

        #region Protected Methods

        protected override string GetStringForOutput()
        {
            StringBuilder sb = new StringBuilder();
            return sb.ToString();
        }

        #endregion Protected Methods
    }
}