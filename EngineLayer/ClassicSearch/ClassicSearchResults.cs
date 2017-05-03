using System.Collections.Generic;
using System.Text;

namespace EngineLayer.ClassicSearch
{
    public class ClassicSearchResults : MetaMorpheusEngineResults
    {

        #region Internal Constructors

        internal ClassicSearchResults(ClassicSearchEngine searchParams) : base(searchParams)
        {
        }

        #endregion Internal Constructors

        #region Public Properties

        public PsmParent[][] OuterPsms { get; set; }

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