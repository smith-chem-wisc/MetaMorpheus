using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Localization
{
    /// <summary>
    /// The probability of each candidate placement of a match's modifications, given its spectrum and the prior learned
    /// from all confident matches, and the localization q-value of its best placement.
    /// </summary>
    public sealed class LocalizationProbabilities
    {
        internal LocalizationProbabilities(IReadOnlyList<LocalizationPlacement> placements, double[] placementProbabilities)
        {
            Placements = placements;
            PlacementProbabilities = placementProbabilities;
            int best = 0;
            for (int c = 1; c < placementProbabilities.Length; c++)
            {
                if (placementProbabilities[c] > placementProbabilities[best])
                {
                    best = c;
                }
            }
            BestPlacementIndex = best;
        }

        /// <summary>The candidate placements, as in the match's <see cref="LocalizationEvidence"/>.</summary>
        public IReadOnlyList<LocalizationPlacement> Placements { get; }

        /// <summary>The posterior probability of each placement; they sum to one.</summary>
        public IReadOnlyList<double> PlacementProbabilities { get; }

        public int BestPlacementIndex { get; }

        public LocalizationPlacement BestPlacement => Placements[BestPlacementIndex];

        /// <summary>The probability that the best placement is the right one.</summary>
        public double BestProbability => PlacementProbabilities[BestPlacementIndex];

        /// <summary>Whether there was more than one placement to choose between.</summary>
        public bool IsLocalizationInQuestion => Placements.Count > 1;

        /// <summary>
        /// The localization q-value of the best placement: the lowest false localization rate (the mean of one minus
        /// the best probability) of any set of confident matches, ranked by best probability, that includes this one.
        /// Null for decoys, for matches that are not confident, and for matches with only one placement.
        /// </summary>
        public double? QValue { get; internal set; }

        /// <summary>
        /// The probability that a modification sits on the residue at <paramref name="residueIndex"/> (for the first and
        /// last residues, on the residue or its peptide terminus): the sum of the probabilities of the placements that
        /// put one there. <paramref name="which"/> limits it to some modifications.
        /// </summary>
        public double SiteProbability(int residueIndex, Func<Modification, bool> which = null)
        {
            double probability = 0;
            for (int c = 0; c < Placements.Count; c++)
            {
                if (Placements[c].Sites.Any(s => s.ResidueIndex == residueIndex && (which == null || which(s.Modification))))
                {
                    probability += PlacementProbabilities[c];
                }
            }
            return probability;
        }
    }
}
