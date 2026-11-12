using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Localization
{
    /// <summary>
    /// The distribution of the relative intensities of experimental peaks that correctly matched fragment ions, as an
    /// empirical cumulative distribution F, learned from the matched peaks of confidently identified unmodified peptides
    /// (An et al., PTMiner, MCP 2019).
    ///
    /// For localization, a peak of relative intensity I counts F(I) toward a placement that explains it and 1 - F(I)
    /// toward one that does not: an intense peak that a placement leaves unexplained is strong evidence against it, a
    /// faint one weak evidence. F is a rank, so any increasing transform of intensity (PTMiner takes the square root)
    /// gives the same F.
    /// </summary>
    public sealed class MatchedPeakIntensityModel
    {
        /// <summary>
        /// Peaks fainter than this, relative to the base peak, are left out of the model and of localization, as in PTMiner.
        /// </summary>
        public const double DefaultMinRelativeIntensity = 0.01;

        private readonly double[] _sorted;

        /// <param name="matchedRelativeIntensities">Relative intensities of correctly matched peaks; those below
        /// <paramref name="minRelativeIntensity"/> are dropped.</param>
        public MatchedPeakIntensityModel(IEnumerable<double> matchedRelativeIntensities, double minRelativeIntensity = DefaultMinRelativeIntensity)
        {
            MinRelativeIntensity = minRelativeIntensity;
            _sorted = matchedRelativeIntensities.Where(i => i >= minRelativeIntensity).OrderBy(i => i).ToArray();
        }

        public double MinRelativeIntensity { get; }

        /// <summary>The number of matched peaks the distribution was learned from.</summary>
        public int TrainingPeakCount => _sorted.Length;

        /// <summary>
        /// F(I), the fraction of training peaks no more intense than <paramref name="relativeIntensity"/>, kept strictly
        /// inside (0, 1) as (count + 0.5) / (n + 1) so that neither F nor 1 - F is ever zero.
        /// </summary>
        public double CumulativeProbability(double relativeIntensity)
        {
            int atOrBelow = UpperBound(relativeIntensity);
            return (atOrBelow + 0.5) / (_sorted.Length + 1);
        }

        /// <summary>The number of training values at or below <paramref name="value"/>.</summary>
        private int UpperBound(double value)
        {
            int lo = 0, hi = _sorted.Length;
            while (lo < hi)
            {
                int mid = (lo + hi) >>> 1;
                if (_sorted[mid] <= value)
                {
                    lo = mid + 1;
                }
                else
                {
                    hi = mid;
                }
            }
            return lo;
        }
    }
}
