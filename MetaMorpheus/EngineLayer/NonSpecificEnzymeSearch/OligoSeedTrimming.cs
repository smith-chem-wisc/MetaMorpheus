using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Fragmentation.Oligo;
using Omics.Modifications;
using Transcriptomics;
using Transcriptomics.Digestion;

namespace EngineLayer.NonSpecificEnzymeSearch
{
    /// <summary>
    /// The nucleic acid half of the non-specific search: which end of a seed oligo is fixed, what the trimmed
    /// candidates weigh, and the candidate itself once one is accepted.
    /// </summary>
    /// <remarks>
    /// <para>The rnase singleN gives, for every position, the longest oligo whose 5' end is fixed there; singleC does
    /// the same for the 3' end. The search decides the other end from the precursor mass, so a seed of length L
    /// stands for its L sub-oligos that keep the fixed end.</para>
    /// <para><b>Masses.</b> A sub-oligo keeping n nucleotides weighs the fixed-end terminus, those n nucleotides,
    /// the modifications on them (and on the fixed terminus), and the terminus the cut leaves at the open end. That
    /// last one is the rnase's remainder (mzLib <see cref="Rnase.ThreePrimeTerminusRemainder"/> for singleN,
    /// <see cref="Rnase.FivePrimeTerminusRemainder"/> for singleC), except for n = L, which keeps the seed's own
    /// terminus. The peptide search gets the same numbers from its b or y ions plus water; for nucleic acids the
    /// fragment caps and termini vary, so they are summed from the nucleotides directly.</para>
    /// <para><b>Modification keys</b> follow <see cref="OligoWithSetMods.AllModsOneIsNterminus"/>: 1 is the 5'
    /// terminus, 2 to L + 1 are nucleotides 1 to L (and the backbone after them), L + 2 is the 3' terminus.</para>
    /// </remarks>
    public static class OligoSeedTrimming
    {
        /// <summary>The end singleN (5') or singleC (3') keeps fixed, or null for any other rnase.</summary>
        public static FragmentationTerminus? FixedTerminus(Rnase rnase) => rnase.CleavageSpecificity switch
        {
            CleavageSpecificity.SingleN => FragmentationTerminus.FivePrime,
            CleavageSpecificity.SingleC => FragmentationTerminus.ThreePrime,
            _ => null,
        };

        /// <summary>The termini a cut can leave at the open end of a seed from this rnase.</summary>
        public static IList<IHasChemicalFormula> OpenEndTermini(Rnase rnase) => FixedTerminus(rnase) == FragmentationTerminus.ThreePrime
            ? rnase.FivePrimeTerminusRemainder
            : rnase.ThreePrimeTerminusRemainder;

        /// <summary>
        /// Fragment types whose mass, less their cap, is exactly the nucleotides they hold plus the fixed terminus, so
        /// that a fragment of n nucleotides points at the sub-oligo of length n. Base-loss ions are excluded: they
        /// lack the base of their last nucleotide. So is M, the intact ion.
        /// </summary>
        public static List<ProductType> TrimmableProductTypes(List<ProductType> dissociationProductTypes, FragmentationTerminus fixedTerminus) =>
            fixedTerminus.GetRnaTerminusSpecificProductTypes()
                .Intersect(dissociationProductTypes)
                .Where(p => !p.IsBaseLoss())
                .ToList();

        /// <summary>
        /// For each trimmable fragment type and open-end terminus, how far the sub-oligo's mass sits above the
        /// fragment's: the open-end terminus less the fragment cap. The non-specific search looks for fragments this
        /// far below the precursor mass.
        /// </summary>
        public static List<double> PrecursorMinusFragmentMasses(IEnumerable<ProductType> trimmableProductTypes, IList<IHasChemicalFormula> openEndTermini) =>
            trimmableProductTypes
                .SelectMany(p => openEndTermini.Select(t => t.MonoisotopicMass - p.GetRnaMassShiftFromProductType()))
                .ToList();

        /// <summary>
        /// The mass of the sub-oligo that keeps <paramref name="length"/> nucleotides at the fixed end, without the
        /// open-end terminus: add one of <see cref="OpenEndTermini"/> to get the precursor mass. Index 0 is unused.
        /// For the full length, add the seed's own open-end terminus instead (or use its MonoisotopicMass).
        /// </summary>
        public static double[] MassesWithoutOpenEndTerminus(OligoWithSetMods seed, FragmentationTerminus fixedTerminus)
        {
            int length = seed.BaseSequence.Length;
            var masses = new double[length + 1];
            bool fivePrimeFixed = fixedTerminus == FragmentationTerminus.FivePrime;
            var mods = seed.AllModsOneIsNterminus;

            double mass = fivePrimeFixed ? seed.FivePrimeTerminus.MonoisotopicMass : seed.ThreePrimeTerminus.MonoisotopicMass;
            if (mods.TryGetValue(fivePrimeFixed ? 1 : length + 2, out Modification fixedTerminalMod))
            {
                mass += fixedTerminalMod.MonoisotopicMass ?? 0;
            }

            for (int n = 1; n <= length; n++)
            {
                int position = fivePrimeFixed ? n : length - n + 1; // one-based position of the nucleotide added
                mass += Nucleotide.GetResidue(seed.BaseSequence[position - 1]).MonoisotopicMass;
                if (mods.TryGetValue(position + 1, out Modification mod))
                {
                    mass += mod.MonoisotopicMass ?? 0;
                }
                masses[n] = mass;
            }
            return masses;
        }

        /// <summary>
        /// The sub-oligo of <paramref name="seed"/> that keeps <paramref name="length"/> nucleotides at the fixed end,
        /// with <paramref name="openEndTerminus"/> at the cut and, if given, <paramref name="openEndMod"/> on it. The
        /// full length keeps the seed's own termini, as the seed's own mass does.
        /// </summary>
        /// <remarks>
        /// Cleavage specificity is None: with no specific rnase to compare the ends with, every candidate is
        /// non-specific, which is also the FDR category it is scored in. As in the peptide search, a trimmed
        /// candidate counts none of its modifications as fixed.
        /// </remarks>
        public static OligoWithSetMods Trim(OligoWithSetMods seed, FragmentationTerminus fixedTerminus, int length,
            IHasChemicalFormula openEndTerminus, Modification openEndMod = null)
        {
            bool fivePrimeFixed = fixedTerminus == FragmentationTerminus.FivePrime;
            int removed = seed.BaseSequence.Length - length;
            int start = fivePrimeFixed ? seed.OneBasedStartResidue : seed.OneBasedStartResidue + removed;
            int end = fivePrimeFixed ? seed.OneBasedEndResidue - removed : seed.OneBasedEndResidue;

            var mods = new Dictionary<int, Modification>();
            foreach (var (key, mod) in seed.AllModsOneIsNterminus)
            {
                // 5' fixed: keep the 5' terminus and nucleotides 1..length (keys 1..length + 1), and the 3' terminus
                // only if nothing was removed. 3' fixed: shift keys down by the nucleotides removed, dropping the 5'
                // terminus unless nothing was removed.
                int newKey = fivePrimeFixed ? key : key - removed;
                bool kept = removed == 0 || (fivePrimeFixed ? key <= length + 1 : newKey >= 2);
                if (kept)
                {
                    mods[newKey] = mod;
                }
            }

            IHasChemicalFormula fivePrimeTerminus = seed.FivePrimeTerminus;
            IHasChemicalFormula threePrimeTerminus = seed.ThreePrimeTerminus;
            if (removed > 0)
            {
                if (fivePrimeFixed)
                {
                    threePrimeTerminus = openEndTerminus;
                }
                else
                {
                    fivePrimeTerminus = openEndTerminus;
                }
            }
            if (openEndMod != null)
            {
                mods[fivePrimeFixed ? length + 2 : 1] = openEndMod;
            }

            return new OligoWithSetMods(seed.NucleicAcid, (RnaDigestionParams)seed.DigestionParams, start, end, 0,
                CleavageSpecificity.None, mods, removed == 0 ? seed.NumFixedMods : 0, fivePrimeTerminus, threePrimeTerminus);
        }
    }
}
