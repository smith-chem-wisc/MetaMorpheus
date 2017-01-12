using InternalLogicEngineLayer;

namespace InternalLogicTaskLayer
{
    public class SearchModeFoSearch
    {
        #region Public Constructors

        public SearchModeFoSearch(SearchMode uu)
        {
            this.sm = uu;
        }

        #endregion Public Constructors

        #region Public Properties

        public SearchMode sm { get; private set; }

        public bool Use { get; set; }

        public string name { get { return sm.FileNameAddition; } }

        #endregion Public Properties
    }
}