using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.FdrAnalysis
{
    /// <summary>
    /// Static math helpers for retention time calibration and z-score computation.
    /// Used by PepAnalysisEngine to convert predicted RT values into a normalized
    /// feature score. Model-agnostic — works for any predictor whose output can be
    /// linearly calibrated against observed RT.
    /// </summary>
    public static class IrtCalibrationHelper
    {
        private const double SentinelZScore = 10.0;
        private const double MinBinStdDev = 0.05;
        private const int BinWidthMinutes = 2;

        /// <summary>
        /// Ordinary least squares regression. Returns (slope, intercept).
        /// Returns (NaN, NaN) if fewer than 2 data points.
        /// </summary>
        public static (double Slope, double Intercept) OlsRegression(
            List<double> x, List<double> y)
        {
            if (x.Count < 2 || x.Count != y.Count)
                return (double.NaN, double.NaN);

            double n = x.Count;
            double sumX = x.Sum();
            double sumY = y.Sum();
            double sumXY = x.Zip(y, (xi, yi) => xi * yi).Sum();
            double sumX2 = x.Sum(xi => xi * xi);

            double denom = n * sumX2 - sumX * sumX;
            if (Math.Abs(denom) < 1e-12)
                return (double.NaN, double.NaN);

            double slope = (n * sumXY - sumX * sumY) / denom;
            double intercept = (sumY - slope * sumX) / n;
            return (slope, intercept);
        }

        /// <summary>
        /// Returns the 2-minute bin index for a given retention time (in minutes).
        /// </summary>
        public static int GetTwoMinuteBin(double retentionTimeMinutes)
            => (int)(retentionTimeMinutes / BinWidthMinutes);

        /// <summary>
        /// Computes per-bin mean and standard deviation of residuals.
        /// Input: list of (observedRT, residual) pairs.
        /// Output: dictionary of bin -> (mean, stddev).
        /// </summary>
        public static Dictionary<int, Tuple<double, double>> ComputeBinStatistics(
            List<(double ObservedRT, double Residual)> residuals)
        {
            var result = new Dictionary<int, Tuple<double, double>>();
            if (!residuals.Any()) return result;

            var grouped = residuals
                .GroupBy(r => GetTwoMinuteBin(r.ObservedRT));

            foreach (var group in grouped)
            {
                var values = group.Select(g => g.Residual).ToList();
                double mean = values.Average();
                double variance = values.Select(v => (v - mean) * (v - mean)).Average();
                double stddev = Math.Sqrt(variance);
                result[group.Key] = Tuple.Create(mean, stddev);
            }

            return result;
        }

        /// <summary>
        /// Raises any bin standard deviation below MinBinStdDev to MinBinStdDev.
        /// Prevents artificially tight bins from inflating z-scores.
        /// Modifies the dictionary in place.
        /// </summary>
        public static void StabilizeBinStDevs(
            Dictionary<int, Tuple<double, double>> binStats)
        {
            foreach (var key in binStats.Keys.ToList())
            {
                if (binStats[key].Item2 < MinBinStdDev)
                    binStats[key] = Tuple.Create(binStats[key].Item1, MinBinStdDev);
            }
        }

        /// <summary>
        /// Computes the z-score of a residual given bin mean and stddev.
        /// Caps at SentinelZScore (10.0). Returns SentinelZScore if stddev is zero.
        /// </summary>
        public static double ComputeZScore(double residual, double mean, double stddev)
        {
            if (stddev <= 0) return SentinelZScore;
            double z = Math.Abs(residual - mean) / stddev;
            return Math.Min(z, SentinelZScore);
        }
    }
}
