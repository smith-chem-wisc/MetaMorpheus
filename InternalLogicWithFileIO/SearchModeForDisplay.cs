using InternalLogicEngineLayer;

namespace InternalLogicTaskLayer
{
    public class SearchModeFoSearch
    {
        public SearchMode sm { get; private set; }

        public bool Use { get; set; }

        public string name { get { return sm.FileNameAddition; } }

        public SearchModeFoSearch(SearchMode uu)
        {
            this.sm = uu;
        }
    }
}