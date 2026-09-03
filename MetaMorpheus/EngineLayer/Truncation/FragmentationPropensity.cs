using System.Collections.Generic;
using MassSpectrometry;
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
    /// First-pass STATIC collisional (HCD/CID) model: enhanced cleavage C-terminal to Asp (strong) and Glu,
    /// and N-terminal to Pro. Weights are relative (an "average" bond = 1.0), so a uniform model reduces
    /// exactly to the standard matched-ion score. Pluggable later for ETD-specific or data-learned
    /// propensities.
    ///
    /// PROVENANCE of the numbers, because this score is what the internal FDR is computed over and the
    /// multipliers therefore set q-values: LeDuc 2014 is cited for the CONCEPT -- weighting a matched ion by
    /// the cleavage propensity of the bond it reports -- and for the three effects being real in collisional
    /// dissociation. The magnitudes 3.0 / 2.0 / 3.0, and the decision to STACK them multiplicatively so that
    /// a D-P bond scores 9.0, are a first pass chosen here, not values read out of that paper. One matched
    /// ion at a D-P bond therefore outweighs nine at average bonds, which is a strong claim for numbers that
    /// have not been fitted. They are the obvious thing to calibrate against a labelled set before this
    /// feature is defaulted on.
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
        /// Whether the static propensity model describes this dissociation type. The three effects it encodes
        /// are collisional: they are properties of a charge-directed backbone cleavage, and they are not
        /// exhibited by the N-Cα cleavage that ECD/ETD produce. Outside these types
        /// <see cref="Score"/> falls back to uniform weights, which is the documented intent -- the model is
        /// described as HCD/CID and "pluggable later for ETD-specific propensities", so applying it to an ETD
        /// scan anyway would push D/E/P weights into an internal FDR that has no basis for them.
        ///
        /// EThcD is deliberately NOT included: it produces b/y and c/z together, so a single scalar weight per
        /// bond cannot be right for both halves of its ladder, and uniform is the honest choice until the
        /// model can be applied per ion series.
        /// </summary>
        public static bool ModelApplies(DissociationType dissociationType) =>
            dissociationType == DissociationType.HCD
            || dissociationType == DissociationType.CID
            || dissociationType == DissociationType.LowCID
            || dissociationType == DissociationType.ISCID;

        /// <summary>
        /// Propensity-weighted match score: sum over matched fragment ions of the cleavage weight at the bond
        /// the ion reports, each scaled by (1 + fractional intensity) — the standard MetaMorpheus per-ion
        /// contribution, multiplied by the cleavage propensity.
        ///
        /// Exactly reproduces <c>MetaMorpheusEngine.CalculatePeptideScore</c> when every weight is 1.0, which
        /// is what licenses comparing an internal-fragment score with a terminal one. Two details are what
        /// make that exact rather than nearly so: a whole-length b/y ion reports no internal bond, so it takes
        /// the average weight of 1.0 rather than being dropped; and diagnostic ions (<see cref="ProductType.D"/>)
        /// count nothing, as they do there.
        /// </summary>
        /// <param name="dissociationType">
        /// The scan's resolved dissociation type. Weights are applied only where
        /// <see cref="ModelApplies"/> says the model holds; everywhere else every bond weighs 1.0 and this
        /// is the standard score.
        /// </param>
        public static double Score(IReadOnlyList<MatchedFragmentIon> matchedIons, string baseSequence,
            double totalIntensity, DissociationType dissociationType)
        {
            bool applyWeights = ModelApplies(dissociationType);
            double score = 0;
            int length = baseSequence.Length;
            foreach (MatchedFragmentIon ion in matchedIons)
            {
                Product product = ion.NeutralTheoreticalProduct;
                if (product.ProductType == ProductType.D)
                {
                    continue; // diagnostic ions count nothing, matching CalculatePeptideScore
                }

                double weight = 1.0;
                if (applyWeights
                    && TryGetCleavageBond(product.Terminus, product.FragmentNumber, length, out int before, out int after))
                {
                    weight = CleavageWeight(baseSequence[before], baseSequence[after]);
                }

                double intensityFraction = totalIntensity > 0 ? ion.Intensity / totalIntensity : 0;
                score += weight * (1 + intensityFraction);
            }
            return score;
        }
    }
}
