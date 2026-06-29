using System.Collections.Generic;
using Omics.Fragmentation;

namespace EngineLayer.Truncation
{
    /// <summary>
    /// Fragmentation-propensity-weighted fragment scoring — the borrowable core of ProSightPD's C-Score
    /// (LeDuc 2014, J. Proteome Res. 13:3231): a matched fragment ion is weighted by the cleavage frequency
    /// of the two residues flanking the backbone bond it reports, so a match at a chemically-favored bond
    /// counts more than a coincidental match at an unfavorable one. This sharpens true-vs-false separation,
    /// which is exactly what the internal-fragment class needs (a real internal fragment is explained by a
    /// chemically-plausible b/y ladder; a random mass-coincidence is not).
    ///
    /// First-pass STATIC HCD/CID model: enhanced cleavage C-terminal to Asp (strong) and Glu, and N-terminal
    /// to Pro. Weights are relative (an "average" bond = 1.0), so with a uniform model this reduces exactly to
    /// the standard matched-ion score. Pluggable later for ETD-specific or data-learned propensities.
    /// </summary>
    public static class FragmentationPropensity
    {
        /// <summary>Relative cleavage weight for the backbone bond between N-side residue <paramref name="before"/>
        /// and C-side residue <paramref name="after"/> (1.0 = average bond).</summary>
        public static double CleavageWeight(char before, char after)
        {
            double weight = 1.0;
            if (before == 'D') weight *= 3.0;        // strong enhanced cleavage C-terminal to aspartate
            else if (before == 'E') weight *= 2.0;   // enhanced cleavage C-terminal to glutamate
            if (after == 'P') weight *= 3.0;         // enhanced cleavage N-terminal to proline
            return weight;
        }

        /// <summary>
        /// Maps a fragment ion (its series + 1-based fragment number) to the 0-based residue indices flanking
        /// the backbone bond it reports, in a peptide of length <paramref name="length"/>. Returns false for a
        /// whole-length ion that reports no internal bond.
        /// </summary>
        public static bool TryGetCleavageBond(FragmentationTerminus terminus, int fragmentNumber, int length,
            out int beforeIndex, out int afterIndex)
        {
            if (terminus == FragmentationTerminus.N)
            {
                beforeIndex = fragmentNumber - 1;   // b/c_k cleaves between residues k and k+1
                afterIndex = fragmentNumber;
            }
            else
            {
                afterIndex = length - fragmentNumber; // y/z_k cleaves between residues L-k and L-k+1
                beforeIndex = afterIndex - 1;
            }
            return beforeIndex >= 0 && afterIndex < length;
        }

        /// <summary>
        /// Propensity-weighted match score: sum over matched fragment ions of the cleavage weight at the bond
        /// the ion reports, each scaled by (1 + fractional intensity) — the standard MetaMorpheus per-ion
        /// contribution, multiplied by the cleavage propensity. A uniform propensity reproduces the standard score.
        /// </summary>
        public static double Score(IReadOnlyList<MatchedFragmentIon> matchedIons, string baseSequence, double totalIntensity)
        {
            double score = 0;
            int length = baseSequence.Length;
            foreach (MatchedFragmentIon ion in matchedIons)
            {
                Product product = ion.NeutralTheoreticalProduct;
                if (!TryGetCleavageBond(product.Terminus, product.FragmentNumber, length, out int before, out int after))
                {
                    continue;
                }
                double weight = CleavageWeight(baseSequence[before], baseSequence[after]);
                double intensityFraction = totalIntensity > 0 ? ion.Intensity / totalIntensity : 0;
                score += weight * (1 + intensityFraction);
            }
            return score;
        }
    }
}
