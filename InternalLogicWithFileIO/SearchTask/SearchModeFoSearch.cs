using InternalLogicEngineLayer;

namespace InternalLogicTaskLayer
{
    public class SearchModeFoSearch
    {
        #region Public Constructors

        public SearchModeFoSearch(SearchMode searchMode)
        {
            SearchMode = searchMode;
        }

        #endregion Public Constructors

        #region Public Properties

        public SearchMode SearchMode { get; private set; }

        public bool Use { get; set; }

        public string Name { get { return SearchMode.FileNameAddition; } }

        #endregion Public Properties
    }
}