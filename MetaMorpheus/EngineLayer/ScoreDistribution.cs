using System;
using System.Threading;

namespace EngineLayer
{
    /// <summary>
    /// Accumulates the scores of every candidate evaluated against one spectrum so that a
    /// percentile of that spectrum's own score distribution can be recovered afterwards.
    /// </summary>
    /// <remarks>
    /// A fixed-width histogram rather than the list of scores it replaces. The search loop is
    /// protein-major and parallel, so scores for a single scan arrive from every thread over the
    /// whole run: a list would need a lock per candidate and would grow without bound, while a
    /// histogram bin can be bumped with <see cref="Interlocked.Increment(ref int)"/> and costs a
    /// fixed 1 kB per scan. The trade is resolution -- a percentile is recovered only to within
    /// half a bin -- and an upper limit, above which scores are clamped into the top bin.
    /// </remarks>
    public class ScoreDistribution
    {
        /// <summary>Default bin width, in score units. The Morpheus score counts matched ions, so a
        /// quarter-ion resolution is finer than the quantity being calibrated.</summary>
        public const double DefaultBinWidth = 0.25;

        /// <summary>Default number of bins, covering scores of 0 to 64 at the default width. A
        /// bottom-up candidate scoring above 64 has matched 64 fragment ions.</summary>
        public const int DefaultBinCount = 256;

        private readonly int[] _bins;
        private int _count;

        public ScoreDistribution(double binWidth = DefaultBinWidth, int binCount = DefaultBinCount)
        {
            if (binWidth <= 0)
                throw new ArgumentOutOfRangeException(nameof(binWidth));
            if (binCount < 1)
                throw new ArgumentOutOfRangeException(nameof(binCount));

            BinWidth = binWidth;
            _bins = new int[binCount];
        }

        public double BinWidth { get; }

        /// <summary>The score at and above which observations are clamped into the top bin.</summary>
        public double ClampThreshold => (_bins.Length - 1) * BinWidth;

        /// <summary>Number of observations recorded.</summary>
        public int Count => Volatile.Read(ref _count);

        /// <summary>
        /// Records one candidate score. Safe to call concurrently from any number of threads.
        /// Negative and non-finite scores are ignored; scores at or above
        /// <see cref="ClampThreshold"/> land in the top bin.
        /// </summary>
        public void Add(double score)
        {
            if (double.IsNaN(score) || score < 0)
                return;

            int bin = double.IsPositiveInfinity(score) ? _bins.Length - 1 : (int)(score / BinWidth);
            if (bin >= _bins.Length)
                bin = _bins.Length - 1;

            Interlocked.Increment(ref _bins[bin]);
            Interlocked.Increment(ref _count);
        }

        /// <summary>
        /// Recovers the requested percentile of the recorded distribution, to within half a bin.
        /// </summary>
        /// <param name="percentile">Fraction in (0, 1]; 0.99 is the 99th percentile.</param>
        /// <param name="value">Centre of the bin holding the observation at that rank.</param>
        /// <returns>
        /// False when the percentile cannot be stated: nothing was recorded, or the rank falls in
        /// the clamped top bin, where the true value is only known to be at or above
        /// <see cref="ClampThreshold"/>. Callers must not substitute a default -- a clamped
        /// distribution yields no answer rather than a wrong one.
        /// </returns>
        public bool TryGetPercentile(double percentile, out double value)
        {
            if (percentile <= 0 || percentile > 1)
                throw new ArgumentOutOfRangeException(nameof(percentile));

            value = double.NaN;
            int total = Count;
            if (total == 0)
                return false;

            // Nearest-rank: the smallest score at or below which at least ceil(p*N) observations fall.
            int targetRank = (int)Math.Ceiling(percentile * total);
            if (targetRank < 1)
                targetRank = 1;

            int cumulative = 0;
            for (int bin = 0; bin < _bins.Length; bin++)
            {
                cumulative += Volatile.Read(ref _bins[bin]);
                if (cumulative < targetRank)
                    continue;

                if (bin == _bins.Length - 1 && ClampThreshold > 0)
                    return false;

                value = (bin + 0.5) * BinWidth;
                return true;
            }

            return false;
        }
    }
}
