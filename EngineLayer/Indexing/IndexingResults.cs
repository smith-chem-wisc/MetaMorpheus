using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.Indexing
{
    public class IndexingResults : MetaMorpheusEngineResults
    {
        public IndexingResults(PeptideWithSetModifications[] peptideIndex, int[][] fragmentIndex, int[][] precursorIndex, IndexingEngine indexParams) : base(indexParams)
        {
            PeptideIndex = peptideIndex;
            FragmentIndex = fragmentIndex;
            PrecursorIndex = precursorIndex;
        }

        public int[][] FragmentIndex { get; private set; }
        public int[][] PrecursorIndex { get; private set; }
        public PeptideWithSetModifications[] PeptideIndex { get; private set; }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("\t\tfragmentIndexDict.Count: " + FragmentIndex.Length);
            if(PrecursorIndex!=null)
            {
                sb.AppendLine("\t\tprecursorIndexDict.Count: " + PrecursorIndex.Length);
            }
            sb.AppendLine("\t\tpeptideIndex.Count: " + PeptideIndex.Length);
            return sb.ToString();
        }
    }
}