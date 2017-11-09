using System.Collections.Generic;
using System.Text;

namespace EngineLayer.Indexing
{
    public class IndexingResults : MetaMorpheusEngineResults
    {
        #region Public Constructors

        public IndexingResults(List<CompactPeptide> peptideIndex, List<int>[] fragmentIndex, IndexingEngine indexParams) : base(indexParams)
        {
            this.PeptideIndex = peptideIndex;
            this.FragmentIndex = fragmentIndex;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<int>[] FragmentIndex { get; private set; }
        public List<CompactPeptide> PeptideIndex { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("\t\tfragmentIndexDict.Count: " + FragmentIndex.Length);
            sb.AppendLine("\t\tpeptideIndex.Count: " + PeptideIndex.Count);
            return sb.ToString();
        }

        #endregion Public Methods
    }
}