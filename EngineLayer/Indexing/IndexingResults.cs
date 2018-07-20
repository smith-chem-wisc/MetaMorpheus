using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.Indexing
{
    public class IndexingResults : MetaMorpheusEngineResults
    {
        public IndexingResults(List<CompactPeptide> peptideIndex, List<int>[] fragmentIndex, IndexingEngine indexParams) : base(indexParams)
        {
            PeptideIndex = peptideIndex;
            FragmentIndex = fragmentIndex;
        }

        public List<int>[] FragmentIndex { get; private set; }
        public List<CompactPeptide> PeptideIndex { get; private set; }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("\t\tfragmentIndexDict.Count: " + FragmentIndex.Length);
            sb.AppendLine("\t\tpeptideIndex.Count: " + PeptideIndex.Count);
            return sb.ToString();
        }
    }
}