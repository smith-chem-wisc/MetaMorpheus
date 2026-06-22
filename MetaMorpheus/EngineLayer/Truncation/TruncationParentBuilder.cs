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
                // A malformed row with no sequence can't become a parent; skip before Split('|') NREs.
                if (string.IsNullOrWhiteSpace(row.FullSequence))
                {
                    continue;
                }

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

        /// <summary>
        /// Returns the combined parent list augmented with a fresh reverse decoy for every target that
        /// does not already have a reverse-decoy counterpart in the input (decision #14). Reuses the
        /// existing mzLib reversal (<see cref="PeptideWithSetModifications.GetReverseDecoyFromTarget"/>),
        /// which locks PTMs to residues; decoys are tagged "DECOY_&lt;accession&gt;". Target parents must
        /// carry a Protein and DigestionParams (the reversal needs them).
        /// </summary>
        public static List<TruncationParent> AddReverseDecoys(IReadOnlyList<TruncationParent> parents)
        {
            var combined = new List<TruncationParent>(parents);

            // Track decoys by (decoy accession + reversed FullSequence) so we add at most one decoy per
            // DISTINCT decoy proteoform — and crucially one per target PROTEOFORM, not one per accession.
            // A protein with several identified proteoforms (e.g. unmodified + oxidized, or pipe-ambiguity
            // alternatives) needs a decoy for each: keying on the accession alone generated a single decoy
            // and left the extra target proteoforms with no decoy counterpart, under-counting decoys and
            // biasing the pooled target-decoy FDR (review finding #6). Seed the set with any decoys already
            // present (e.g. decoy hits inherited from Pass 1) so those are not regenerated.
            var existingDecoys = parents
                .Where(p => p.IsDecoy)
                .Select(p => DecoyKey(p.ProteinAccession, p.Proteoform.FullSequence))
                .ToHashSet();

            foreach (TruncationParent target in parents.Where(p => !p.IsDecoy))
            {
                string decoyAccession = "DECOY_" + target.ProteinAccession;
                PeptideWithSetModifications decoyProteoform =
                    target.Proteoform.GetReverseDecoyFromTarget(new int[target.Proteoform.BaseSequence.Length]);

                if (!existingDecoys.Add(DecoyKey(decoyAccession, decoyProteoform.FullSequence)))
                {
                    continue; // an identical decoy proteoform already exists (Pass 1, or a duplicate target)
                }

                combined.Add(new TruncationParent(decoyProteoform, decoyAccession, target.OriginatingId, isDecoy: true));
            }

            return combined;
        }

        // Identity of a decoy proteoform for de-duplication: its accession plus its (reversed) full sequence.
        private static string DecoyKey(string accession, string fullSequence) => accession + "|" + fullSequence;
    }
}
