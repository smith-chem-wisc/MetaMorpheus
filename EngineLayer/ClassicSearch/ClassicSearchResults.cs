using System.Text;

namespace EngineLayer.ClassicSearch
{
    public class ClassicSearchResults : MyResults
    {

        #region Internal Constructors

        internal ClassicSearchResults(ClassicSearchEngine searchParams) : base(searchParams)
        {
        }

        #endregion Internal Constructors

        #region Public Properties

        public PsmClassic[][] OuterPsms { get; set; }

        #endregion Public Properties

        #region Protected Properties

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            return sb.ToString();
        }

        #endregion Protected Properties

    }
}