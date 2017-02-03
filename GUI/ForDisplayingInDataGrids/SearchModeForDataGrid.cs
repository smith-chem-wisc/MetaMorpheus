using EngineLayer;

namespace MetaMorpheusGUI
{
    class SearchModeForDataGrid
    {

        #region Public Fields

        public readonly SearchMode searchMode;

        #endregion Public Fields

        #region Public Constructors

        public SearchModeForDataGrid(SearchMode searchMode)
        {
            this.searchMode = searchMode;
        }

        #endregion Public Constructors

        #region Public Properties

        public bool Use { get; set; }

        public string Name { get { return searchMode.FileNameAddition; } }

        #endregion Public Properties

    }
}