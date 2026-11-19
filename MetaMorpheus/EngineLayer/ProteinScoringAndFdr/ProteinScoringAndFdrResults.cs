using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public class ProteinScoringAndFdrResults : MetaMorpheusEngineResults
    {
        public List<ProteinGroup> SortedAndScoredProteinGroups;

        /// <summary>The threshold of the search's identification filter, which the protein count is held to.</summary>
        public double QValueThreshold { get; init; } = 0.01;

        public ProteinScoringAndFdrResults(ProteinScoringAndFdrEngine proteinAnalysisEngine) : base(proteinAnalysisEngine)
        {
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.Append($"Number of proteins with q-value <= {QValueThreshold.ToString(CultureInfo.InvariantCulture)}: " +
                SortedAndScoredProteinGroups.Count(b => b.QValue <= QValueThreshold));
            return sb.ToString();
        }
    }
}
