using System.Collections.Generic;
using System.Linq;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer.Truncation
{
    /// <summary>Outcome of a successful chop: the deduced truncated form and how it was reached.</summary>
    public class ChopResult
    {
        public ChopResult(PeptideWithSetModifications truncatedForm, int residuesChopped, FragmentationTerminus terminusChopped, int notch)
        {
            TruncatedForm = truncatedForm;
            ResiduesChopped = residuesChopped;
            TerminusChopped = terminusChopped;
            Notch = notch;
        }

        public PeptideWithSetModifications TruncatedForm { get; }
        public int ResiduesChopped { get; }
        public FragmentationTerminus TerminusChopped { get; }
        public int Notch { get; }
    }

    /// <summary>
    /// Pass 3 chopping algorithm (docs/Truncation-Search.md decisions #9, #10). Removes residues one at a time
    /// from the indicated terminus — each step also dropping any PTM locked to that residue and the
    /// terminal mod on it (so N-terminal acetyl-Ala leaves as one step) — until the deduced truncated
    /// mass matches the observed precursor mass at an allowed notch of the supplied acceptor, or no
    /// residues remain. Only normal backbone bonds are cut; there is no mass-shift absorption (#9).
    /// </summary>
    public static class ProteoformChopper
    {
        /// <summary>
        /// Chops residues from <paramref name="terminusToChop"/> of <paramref name="parent"/> until the
        /// truncated monoisotopic mass matches <paramref name="targetMass"/> within an allowed notch of
        /// <paramref name="massDiffAcceptor"/> (which carries the precursor tolerance). Returns the first
        /// match, or null if no integer-residue chop matches. At least one residue is always removed.
        /// Chopping is monotonic in mass, so the loop stops early once the truncated form drops below the
        /// lowest mass the acceptor could accept — no point building the remaining forms.
        /// </summary>
        /// <param name="fixedModifications">The search's fixed mods, used to keep the truncated form's fixed/variable
        /// mod split right. When null, the split is inferred from the parent's own counts.</param>
        public static ChopResult ChopUntilMassMatches(PeptideWithSetModifications parent, FragmentationTerminus terminusToChop,
            double targetMass, MassDiffAcceptor massDiffAcceptor, IReadOnlyCollection<Modification> fixedModifications = null)
        {
            int length = parent.BaseSequence.Length;
            double lowestAcceptableMass = LowestAcceptableMass(targetMass, massDiffAcceptor);

            for (int chopped = 1; chopped <= length - 1; chopped++)
            {
                int chopFromN = terminusToChop == FragmentationTerminus.N ? chopped : 0;
                int chopFromC = terminusToChop == FragmentationTerminus.C ? chopped : 0;

                PeptideWithSetModifications truncated = BuildTruncated(parent, chopFromN, chopFromC, fixedModifications);
                double mass = truncated.MonoisotopicMass;

                int notch = massDiffAcceptor.Accepts(targetMass, mass);
                if (notch >= 0)
                {
                    return new ChopResult(truncated, chopped, terminusToChop, notch);
                }

                if (mass < lowestAcceptableMass)
                {
                    // Already lighter than anything the acceptor could take, and from the first chop onwards
                    // the mass only falls (the terminal mod is gone after step 1, so every further step just
                    // removes a residue). No deeper chop can match; stop instead of building the rest.
                    break;
                }
            }

            return null;
        }

        /// <summary>
        /// Lowest theoretical mass <paramref name="massDiffAcceptor"/> could accept for an observed
        /// <paramref name="targetMass"/> — the acceptor's own lower bound across all its notches, which is
        /// what bounds the chopping loops. Returns <see cref="double.NegativeInfinity"/> when the acceptor
        /// reports no bounded interval (open searches), which simply disables the early exit.
        /// </summary>
        private static double LowestAcceptableMass(double targetMass, MassDiffAcceptor massDiffAcceptor)
        {
            double lowest = double.PositiveInfinity;

            foreach (AllowedIntervalWithNotch interval in massDiffAcceptor.GetAllowedPrecursorMassIntervalsFromObservedMass(targetMass))
            {
                if (interval.Minimum < lowest)
                {
                    lowest = interval.Minimum;
                }
            }

            return double.IsPositiveInfinity(lowest) ? double.NegativeInfinity : lowest;
        }

        /// <summary>
        /// Builds the proteoform obtained by removing <paramref name="chopFromN"/> residues from the
        /// N-terminus and <paramref name="chopFromC"/> from the C-terminus of <paramref name="parent"/>,
        /// re-indexing the surviving modifications onto the shorter form. A terminal mod survives only if
        /// no residue was chopped from its terminus (so chopping a terminal residue removes its terminal
        /// mod too, #9/#10). Coordinates are reported within the parent protein (#13).
        /// </summary>
        private static PeptideWithSetModifications BuildTruncated(PeptideWithSetModifications parent, int chopFromN, int chopFromC,
            IReadOnlyCollection<Modification> fixedModifications)
        {
            int parentLength = parent.BaseSequence.Length;
            int newLength = parentLength - chopFromN - chopFromC;
            int newStart = parent.OneBasedStartResidueInProtein + chopFromN;
            int newEnd = parent.OneBasedEndResidueInProtein - chopFromC;

            Dictionary<int, Modification> parentMods = parent.AllModsOneIsNterminus;
            var newMods = new Dictionary<int, Modification>();

            // N-terminal mod survives only if the N-terminus was not chopped.
            if (chopFromN == 0 && parentMods.TryGetValue(1, out Modification nTermMod))
            {
                newMods[1] = nTermMod;
            }

            // Residue mods: new residue r corresponds to parent residue (r + chopFromN); mod keys are residue + 1.
            for (int r = 1; r <= newLength; r++)
            {
                int parentKey = (r + chopFromN) + 1;
                if (parentMods.TryGetValue(parentKey, out Modification residueMod))
                {
                    newMods[r + 1] = residueMod;
                }
            }

            // C-terminal mod survives only if the C-terminus was not chopped.
            if (chopFromC == 0 && parentMods.TryGetValue(parentLength + 2, out Modification cTermMod))
            {
                newMods[newLength + 2] = cTermMod;
            }

            // The peptide Description is emitted verbatim in the standard psmtsv "Description" column
            // (decision #13), mirroring existing proteolysis-product descriptors like "chain(2-121)".
            // A clean initiator-Met excision is labeled as such, not as a 1-residue N-terminal truncation.
            string truncationDescription = ClassifyChop(parent, chopFromN, chopFromC) + $"({newStart}-{newEnd})";

            return new PeptideWithSetModifications(parent.Protein, parent.DigestionParams, newStart, newEnd,
                CleavageSpecificity.Full, truncationDescription, 0, newMods, CountFixedMods(parent, newMods, fixedModifications));
        }

        /// <summary>
        /// How many of the surviving mods are fixed, so NumVariableMods (= NumMods - NumFixedMods) stays right on the
        /// truncated form. With the search's fixed mods, count the survivors among them; without, a parent with no
        /// fixed mods (usual in top-down) or only fixed mods is exact, and anything else is capped by the parent's count.
        /// </summary>
        private static int CountFixedMods(PeptideWithSetModifications parent, Dictionary<int, Modification> newMods,
            IReadOnlyCollection<Modification> fixedModifications)
        {
            if (fixedModifications != null)
            {
                return newMods.Values.Count(fixedModifications.Contains);
            }

            if (parent.NumFixedMods == 0)
            {
                return 0;
            }

            return parent.NumFixedMods == parent.AllModsOneIsNterminus.Count
                ? newMods.Count
                : System.Math.Min(parent.NumFixedMods, newMods.Count);
        }

        /// <summary>Builds the internal/terminal truncated form for a specific (chopFromN, chopFromC) — used by
        /// the internal search, which finds mass-matching spans by prefix-sum arithmetic and only then builds
        /// the (few) matches rather than constructing every span.</summary>
        public static PeptideWithSetModifications BuildTruncatedForm(PeptideWithSetModifications parent, int chopFromN, int chopFromC)
            => BuildTruncated(parent, chopFromN, chopFromC, null);

        /// <summary>
        /// Enumerates INTERNAL truncations of <paramref name="parent"/> (both termini lost: chopFromN≥1 AND
        /// chopFromC≥1, keeping ≥1 residue) whose mass matches <paramref name="targetMass"/> at an allowed
        /// notch. Keeps the full-length proteoform as the parent (no DB-wide subsequence search). For a fixed
        /// N-chop the mass decreases monotonically in the C-chop, so the inner loop stops once the form is
        /// lighter than any acceptable target — near-linear in length. The caller scores each candidate.
        /// </summary>
        public static List<ChopResult> ChopInternalCandidates(PeptideWithSetModifications parent, double targetMass,
            MassDiffAcceptor massDiffAcceptor)
        {
            int length = parent.BaseSequence.Length;
            var results = new List<ChopResult>();
            double lowestAcceptableMass = LowestAcceptableMass(targetMass, massDiffAcceptor);

            for (int chopFromN = 1; chopFromN <= length - 2; chopFromN++)
            {
                for (int chopFromC = 1; chopFromN + chopFromC <= length - 1; chopFromC++)
                {
                    PeptideWithSetModifications truncated = BuildTruncated(parent, chopFromN, chopFromC, null);
                    double mass = truncated.MonoisotopicMass;

                    int notch = massDiffAcceptor.Accepts(targetMass, mass);
                    if (notch >= 0)
                    {
                        results.Add(new ChopResult(truncated, chopFromN + chopFromC, FragmentationTerminus.None, notch));
                    }
                    else if (mass < lowestAcceptableMass)
                    {
                        break; // already too light; removing more C-terminal residues cannot reach the target
                    }
                }
            }

            return results;
        }

        /// <summary>
        /// Classifies a chop into a Description label (#13). A chop from both termini is an internal truncation.
        /// A C-terminal chop is a C-terminal truncation.
        /// An N-terminal chop is a true N-terminal truncation UNLESS it is a clean removal of exactly the
        /// initiator methionine (one residue, from a parent that starts at protein position 1 with Met),
        /// which is the canonical NME form — reported as "N-terminal Met excision", so NME is not miscounted as
        /// a truncation. There is no "NME + acetylation" label: chopping drops the N-terminal mod (key 1), and
        /// the chopper adds no mass, so an acetyl on the new N-terminus can never be produced here.
        /// </summary>
        private static string ClassifyChop(PeptideWithSetModifications parent, int chopFromN, int chopFromC)
        {
            if (chopFromN > 0 && chopFromC > 0)
            {
                return TruncationPass3.InternalTruncation;
            }
            if (chopFromN == 0)
            {
                return TruncationPass3.CTerminalTruncation;
            }

            bool isInitiatorMetExcision = chopFromN == 1 && chopFromC == 0
                && parent.OneBasedStartResidueInProtein == 1
                && parent.BaseSequence.Length > 0 && parent.BaseSequence[0] == 'M';
            return isInitiatorMetExcision ? TruncationPass3.NTerminalMetExcision : TruncationPass3.NTerminalTruncation;
        }
    }
}
