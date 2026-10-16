using System;

namespace EngineLayer
{
    /// <summary>
    /// Per-spectrum score calibration after Tailor (Yilmaz et al., J. Proteome Res. 2020): a match is
    /// divided by a high percentile of the score distribution of all candidates evaluated against the
    /// same spectrum, so that scores from spectra of differing quality become comparable without any
    /// trained model.
    /// </summary>
    public static class TailorScoreCalculator
    {
        /// <summary>The percentile of the per-spectrum distribution used as the divisor.</summary>
        public const double CalibrationPercentile = 0.99;

        /// <summary>
        /// Fewest candidates a spectrum must have before its percentile is used.
        /// </summary>
        /// <remarks>
        /// The divisor is the ceil(0.99*N)-th smallest of N candidate scores, so it approaches the
        /// best score as N falls: at N=100 it is the second largest and the ratio is pinned near 1.
        /// A spectrum below this count is reported as uncalibrated rather than as a ratio that is
        /// arithmetically valid and carries no information.
        /// </remarks>
        public const int MinimumCandidates = 500;

        /// <summary>
        /// Calibrates <paramref name="score"/> against the distribution its spectrum produced.
        /// </summary>
        /// <returns>
        /// The calibrated score, or NaN when the spectrum cannot be calibrated -- too few
        /// candidates, a distribution whose percentile is unresolved, or a non-positive divisor.
        /// NaN is written as a blank cell; it is never silently replaced by the raw score, which
        /// would make an uncalibrated row indistinguishable from a calibrated one.
        /// </returns>
        public static double Calculate(double score, ScoreDistribution distribution,
            int minimumCandidates = MinimumCandidates)
        {
            if (distribution == null || distribution.Count < minimumCandidates)
                return double.NaN;

            if (!distribution.TryGetPercentile(CalibrationPercentile, out double divisor))
                return double.NaN;

            if (divisor <= 0 || !double.IsFinite(divisor))
                return double.NaN;

            return score / divisor;
        }
    }
}
