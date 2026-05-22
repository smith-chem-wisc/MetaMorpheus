using System.Collections.Generic;
using System.Linq;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer.Truncation
{
    /// <summary>
    /// Minimal representation of a deduped Pass 1 proteoform-level row (the AllProteoforms.psmtsv set,
    /// decision #1) used to seed Pass 2 parents. <see cref="FullSequence"/> may carry pipe-separated
    /// alternatives, which <see cref="TruncationParentBuilder"/> expands into separate parents (#2).
    /// </summary>
    public class Pass1ProteoformRow
    {
        public string FullSequence { get; init; }
        public string ProteinAccession { get; init; }

        /// <summary>PEP q-value, or null when PEP was not computed (too few proteoforms to train).</summary>
        public double? PepQValue { get; init; }

        public double NotchQValue { get; init; }
        public bool IsDecoy { get; init; }
    }

    /// <summary>
    /// Turns deduped Pass 1 proteoform rows into the combined target+decoy parent list for Pass 2
    /// (01_Architecture.md decisions #2, #3). Applies the permissive parent inclusion threshold and
    /// expands pipe-ambiguous FullSequences, tagging every alternative with its originating row so
    /// downstream counting does not double-count.
    /// </summary>
    public static class TruncationParentBuilder
    {
        public static List<TruncationParent> BuildParents(
            IEnumerable<Pass1ProteoformRow> rows,
            double parentQValueThreshold,
            Dictionary<string, Modification> knownModifications = null)
        {
            var parents = new List<TruncationParent>();

            foreach (Pass1ProteoformRow row in rows)
            {
                if (!PassesParentFilter(row, parentQValueThreshold))
                {
                    continue;
                }

                // Pipe-ambiguous: each alternative is a separate parent tied back to the same row (#2).
                foreach (string alternative in row.FullSequence.Split('|'))
                {
                    var proteoform = new PeptideWithSetModifications(alternative, knownModifications);
                    parents.Add(new TruncationParent(proteoform, row.ProteinAccession, row, row.IsDecoy));
                }
            }

            return parents;
        }

        /// <summary>
        /// Permissive parent inclusion (#3): PEP q-value ≤ threshold when PEP was computed; otherwise
        /// notch q-value ≤ threshold. The same filter applies to targets and decoys.
        /// </summary>
        public static bool PassesParentFilter(Pass1ProteoformRow row, double parentQValueThreshold) =>
            row.PepQValue.HasValue
                ? row.PepQValue.Value <= parentQValueThreshold
                : row.NotchQValue <= parentQValueThreshold;
    }
}
