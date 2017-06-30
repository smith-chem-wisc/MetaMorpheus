using EngineLayer;

namespace MetaMorpheusGUI
{
    class SearchModeForDataGrid
    {

        #region Public Fields

        public readonly MassDiffAcceptor searchMode;

        #endregion Public Fields

        #region Public Constructors

        public SearchModeForDataGrid(MassDiffAcceptor searchMode)
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