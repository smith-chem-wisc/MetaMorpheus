using System.Collections.Generic;

namespace InternalLogic
{
    public class IndexResults : MyResults
    {
        public IndexResults(List<CompactPeptide> myDictionary, Dictionary<float, List<int>> myFragmentDictionary, IndexEngine indexParams) : base(indexParams)
        {
            this.peptideIndex = myDictionary;
            this.fragmentIndexDict = myFragmentDictionary;
        }

        public Dictionary<float, List<int>> fragmentIndexDict { get; private set; }
        public List<CompactPeptide> peptideIndex { get; private set; }
    }
}