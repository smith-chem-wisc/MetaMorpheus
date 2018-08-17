using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.Indexing
{
    public class IndexingResultsWithNGlyco : MetaMorpheusEngineResults
    {
        public IndexingResultsWithNGlyco(List<CompactPeptide> peptideIndex, List<int>[] fragmentIndex, List<int>[] fragmentIndexNGly,IndexingEngine indexParams) : base(indexParams)
        {
            PeptideIndex = peptideIndex;
            FragmentIndex = fragmentIndex;
            FragmentIndexNgly = fragmentIndexNGly;
        }

        public List<int>[] FragmentIndex { get; private set; }
        public List<CompactPeptide> PeptideIndex { get; private set; }
        public List<int>[] FragmentIndexNgly { get; private set; }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("\t\tfragmentIndexDict.Count: " + FragmentIndex.Length);
            sb.AppendLine("\t\tpeptideIndex.Count: " + PeptideIndex.Count);
            sb.AppendLine("\t\tfrgmentIndexNGlyDict.Count: " + FragmentIndexNgly.Length);
            return sb.ToString();
        }
    }
}
