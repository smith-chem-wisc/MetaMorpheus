using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.Indexing
{
    public class IndexingResults : MetaMorpheusEngineResults
    {
        public IndexingResults(List<PeptideWithSetModifications> peptideIndex, List<int>[] fragmentIndex, List<int>[] precursorIndex, IndexingEngine indexParams) : base(indexParams)
        {
            PeptideIndex = peptideIndex;
            FragmentIndex = fragmentIndex;
            PrecursorIndex = precursorIndex;
        }

        public List<int>[] FragmentIndex { get; private set; }
        public List<int>[] PrecursorIndex { get; private set; }
        public List<PeptideWithSetModifications> PeptideIndex { get; private set; }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("\t\tfragmentIndexDict.Count: " + FragmentIndex.Length);
            if(PrecursorIndex!=null)
            {
                sb.AppendLine("\t\tprecursorIndexDict.Count: " + PrecursorIndex.Length);
            }
            sb.AppendLine("\t\tpeptideIndex.Count: " + PeptideIndex.Count);
            return sb.ToString();
        }
    }
}