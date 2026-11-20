using MathNet.Numerics.Distributions;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Localization
{
    /// <summary>
    /// How <see cref="MatchedPeakIntensityModel"/> represents the distribution of matched-peak intensities.
    /// </summary>
    public enum MatchedPeakIntensityDistribution
    {
        /// <summary>The empirical cumulative distribution of the training intensities; assumes no shape.</summary>
        Empirical,

        /// <summary>
        /// A lognormal distribution fitted to the training intensities, as PTMiner does (it fits the square root of
        /// relative intensity, which is lognormal exactly when intensity is, with the same cumulative distribution).
        /// </summary>
        LogNormal
    }

    /// <summary>
    /// The distribution of the relative intensities of experimental peaks that correctly matched fragment ions, as a
    /// cumulative distribution F, learned from the matched peaks of confidently identified unmodified peptides
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

        /// <summary>F is kept at least this far from 0 and 1, so that neither log F nor log(1 - F) is ever infinite.</summary>
        private const double LogNormalClamp = 1e-6;

        private readonly double[] _sorted;
        private readonly double _logMean;
        private readonly double _logStandardDeviation;

        /// <param name="matchedRelativeIntensities">Relative intensities of correctly matched peaks; those below
        /// <paramref name="minRelativeIntensity"/> are dropped.</param>
        public MatchedPeakIntensityModel(IEnumerable<double> matchedRelativeIntensities, double minRelativeIntensity = DefaultMinRelativeIntensity,
            MatchedPeakIntensityDistribution distribution = MatchedPeakIntensityDistribution.Empirical)
        {
            MinRelativeIntensity = minRelativeIntensity;
            Distribution = distribution;
            _sorted = matchedRelativeIntensities.Where(i => i >= minRelativeIntensity).OrderBy(i => i).ToArray();

            if (distribution == MatchedPeakIntensityDistribution.LogNormal && _sorted.Length > 1)
            {
                double[] logs = _sorted.Select(i => Math.Log(i)).ToArray();
                _logMean = logs.Average();
                _logStandardDeviation = Math.Sqrt(logs.Sum(l => (l - _logMean) * (l - _logMean)) / (logs.Length - 1));
            }
        }

        public double MinRelativeIntensity { get; }

        public MatchedPeakIntensityDistribution Distribution { get; }

        /// <summary>The number of matched peaks the distribution was learned from.</summary>
        public int TrainingPeakCount => _sorted.Length;

        /// <summary>
        /// F(I), the probability that a correctly matched peak is no more intense than <paramref name="relativeIntensity"/>,
        /// kept strictly inside (0, 1) so that neither F nor 1 - F is ever zero. Empirically it is the fraction of training
        /// peaks at or below I, as (count + 0.5) / (n + 1).
        /// </summary>
        public double CumulativeProbability(double relativeIntensity)
        {
            if (Distribution == MatchedPeakIntensityDistribution.LogNormal)
            {
                if (_logStandardDeviation <= 0 || relativeIntensity <= 0)
                {
                    return relativeIntensity <= 0 ? LogNormalClamp : 0.5;
                }
                double f = Normal.CDF(_logMean, _logStandardDeviation, Math.Log(relativeIntensity));
                return Math.Clamp(f, LogNormalClamp, 1 - LogNormalClamp);
            }

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
