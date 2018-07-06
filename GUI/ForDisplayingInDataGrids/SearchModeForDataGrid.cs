using EngineLayer;

namespace MetaMorpheusGUI
{
    internal class SearchModeForDataGrid
    {
        public readonly MassDiffAcceptor SearchMode;

        public SearchModeForDataGrid(MassDiffAcceptor searchMode)
        {
            SearchMode = searchMode;
        }

        public bool Use { get; set; }

        public string Name { get { return SearchMode.ToString(); } }
    }
}