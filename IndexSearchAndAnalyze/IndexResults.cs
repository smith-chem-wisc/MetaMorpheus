using System.Collections.Generic;

namespace IndexSearchAndAnalyze
{
    public class IndexResults : MyResults
    {
        public IndexResults(List<CompactPeptide> myDictionary, Dictionary<float, List<int>> myFragmentDictionary)
        {
            this.peptideIndex = myDictionary;
            this.fragmentIndexDict = myFragmentDictionary;
        }

        public Dictionary<float, List<int>> fragmentIndexDict { get; private set; }
        public List<CompactPeptide> peptideIndex { get; private set; }
    }
}
