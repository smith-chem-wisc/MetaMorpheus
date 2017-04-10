using System.Collections.Generic;
using System.Text;

namespace EngineLayer.Indexing
{
    public class IndexingResults : MetaMorpheusEngineResults
    {

        #region Public Constructors

        public IndexingResults(List<CompactPeptide> peptideIndex, Dictionary<float, List<int>> fragmentIndexDict, IndexingEngine indexParams) : base(indexParams)
        {
            this.PeptideIndex = peptideIndex;
            this.FragmentIndexDict = fragmentIndexDict;
        }

        #endregion Public Constructors

        #region Public Properties

        public Dictionary<float, List<int>> FragmentIndexDict { get; private set; }
        public List<CompactPeptide> PeptideIndex { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append(base.ToString());
            sb.AppendLine("\t\tfragmentIndexDict.Count: " + FragmentIndexDict.Count);
            sb.AppendLine("\t\tpeptideIndex.Count: " + PeptideIndex.Count);
            return sb.ToString();
        }

        #endregion Public Methods

    }
}