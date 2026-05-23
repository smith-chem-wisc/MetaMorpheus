using System.Collections.Generic;
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
    /// Pass 3 chopping algorithm (01_Architecture.md decisions #9, #10). Removes residues one at a time
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
        /// </summary>
        public static ChopResult ChopUntilMassMatches(PeptideWithSetModifications parent, FragmentationTerminus terminusToChop,
            double targetMass, MassDiffAcceptor massDiffAcceptor)
        {
            int length = parent.BaseSequence.Length;

            for (int chopped = 1; chopped <= length - 1; chopped++)
            {
                int chopFromN = terminusToChop == FragmentationTerminus.N ? chopped : 0;
                int chopFromC = terminusToChop == FragmentationTerminus.C ? chopped : 0;

                PeptideWithSetModifications truncated = BuildTruncated(parent, chopFromN, chopFromC);

                int notch = massDiffAcceptor.Accepts(targetMass, truncated.MonoisotopicMass);
                if (notch >= 0)
                {
                    return new ChopResult(truncated, chopped, terminusToChop, notch);
                }
            }

            return null;
        }

        /// <summary>
        /// Builds the proteoform obtained by removing <paramref name="chopFromN"/> residues from the
        /// N-terminus and <paramref name="chopFromC"/> from the C-terminus of <paramref name="parent"/>,
        /// re-indexing the surviving modifications onto the shorter form. A terminal mod survives only if
        /// no residue was chopped from its terminus (so chopping a terminal residue removes its terminal
        /// mod too, #9/#10). Coordinates are reported within the parent protein (#13).
        /// </summary>
        private static PeptideWithSetModifications BuildTruncated(PeptideWithSetModifications parent, int chopFromN, int chopFromC)
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
            string truncationDescription = ClassifyChop(parent, chopFromN, chopFromC, newMods) + $"({newStart}-{newEnd})";

            return new PeptideWithSetModifications(parent.Protein, parent.DigestionParams, newStart, newEnd,
                CleavageSpecificity.Full, truncationDescription, 0, newMods, newMods.Count);
        }

        /// <summary>
        /// Classifies a chop into a Description label (#13). A C-terminal chop is a C-terminal truncation.
        /// An N-terminal chop is a true N-terminal truncation UNLESS it is a clean removal of exactly the
        /// initiator methionine (one residue, from a parent that starts at protein position 1 with Met),
        /// which is the canonical NME form — reported as "N-terminal Met excision" (or "... + acetylation"
        /// when the new N-terminus carries an acetyl), so NME is not miscounted as a truncation.
        /// </summary>
        private static string ClassifyChop(PeptideWithSetModifications parent, int chopFromN, int chopFromC,
            Dictionary<int, Modification> newMods)
        {
            if (chopFromN == 0)
            {
                return TruncationPass3.CTerminalTruncation;
            }

            bool isInitiatorMetExcision = chopFromN == 1 && chopFromC == 0
                && parent.OneBasedStartResidueInProtein == 1
                && parent.BaseSequence.Length > 0 && parent.BaseSequence[0] == 'M';
            if (isInitiatorMetExcision)
            {
                return HasNTerminalAcetyl(newMods)
                    ? TruncationPass3.NTerminalMetExcisionPlusAcetyl
                    : TruncationPass3.NTerminalMetExcision;
            }

            return TruncationPass3.NTerminalTruncation;
        }

        /// <summary>True if the truncated form carries an acetyl at its N-terminus (mod key 1) or first
        /// residue (mod key 2) — used to label NME + acetylation (#13).</summary>
        private static bool HasNTerminalAcetyl(Dictionary<int, Modification> newMods)
        {
            foreach (KeyValuePair<int, Modification> kv in newMods)
            {
                if (kv.Key <= 2 && kv.Value != null
                    && (kv.Value.OriginalId ?? string.Empty).IndexOf("acetyl", System.StringComparison.OrdinalIgnoreCase) >= 0)
                {
                    return true;
                }
            }
            return false;
        }
    }
}
