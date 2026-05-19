using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.FdrAnalysis
{
    public static class IrtCalibrationHelper
    {
        public const double MaxZScore = 10.0;
        private const int MinCalibrationPeptides = 100;

        /// <summary>
        /// Ordinary least squares linear regression.
        /// Returns (double.NaN, double.NaN) if fewer than MinCalibrationPeptides points.
        /// Returns (slope, intercept) where y = slope * x + intercept.
        /// </summary>
        public static (double Slope, double Intercept) LinearRegression(
            IReadOnlyList<double> x, IReadOnlyList<double> y)
        {
            int n = x.Count;
            if (n < MinCalibrationPeptides)
                return (double.NaN, double.NaN);

            double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
            for (int i = 0; i < n; i++)
            {
                sumX += x[i];
                sumY += y[i];
                sumXY += x[i] * y[i];
                sumX2 += x[i] * x[i];
            }

            double denom = n * sumX2 - sumX * sumX;
            if (Math.Abs(denom) < 1e-15)
                return (double.NaN, double.NaN);

            double slope = (n * sumXY - sumX * sumY) / denom;
            double intercept = (sumY - slope * sumX) / n;
            return (slope, intercept);
        }

        /// <summary>
        /// Groups residuals by 2-minute RT bin and computes (mean, stddev) per bin.
        /// Bin key = (int)(2 * Math.Round(observedRT / 2d, 0))
        /// </summary>
        public static Dictionary<int, Tuple<double, double>> ComputeBinStatistics(
            Dictionary<int, List<double>> residualsByBin)
        {
            var result = new Dictionary<int, Tuple<double, double>>();
            foreach (var kvp in residualsByBin)
            {
                double mean = kvp.Value.Mean();
                double stddev = kvp.Value.StandardDeviation();
                result[kvp.Key] = Tuple.Create(mean, stddev);
            }
            return result;
        }

        /// <summary>
        /// Replaces per-bin stddevs that are too tight (&lt; 0.05) or too wide (ratio to global &gt; 3)
        /// with the global stddev. Modifies averagesCommaStandardDeviations in place.
        /// This mirrors the stabilization logic in PepAnalysisEngine for hydrophobicity.
        /// </summary>
        public static void StabilizeBinStDevs(
            Dictionary<int, Tuple<double, double>> averagesCommaStandardDeviations,
            double globalStdDev)
        {
            foreach (var key in averagesCommaStandardDeviations.Keys.ToList())
            {
                var tuple = averagesCommaStandardDeviations[key];
                double stddev = tuple.Item2;
                if (stddev < 0.05 || (stddev / globalStdDev) > 3)
                {
                    averagesCommaStandardDeviations[key] = Tuple.Create(tuple.Item1, globalStdDev);
                }
            }
        }

        /// <summary>
        /// Computes |residual - binMean| / binStddev, clamped to [0, MaxZScore].
        /// Returns MaxZScore for NaN, infinity, or values exceeding MaxZScore.
        /// </summary>
        public static double ComputeZScore(double residual, double binMean, double binStdDev)
        {
            if (binStdDev == 0 || double.IsNaN(binStdDev) || double.IsInfinity(binStdDev))
                return MaxZScore;

            double z = Math.Abs(residual - binMean) / binStdDev;

            if (double.IsNaN(z) || double.IsInfinity(z) || z > MaxZScore)
                return MaxZScore;

            return z;
        }

        /// <summary>
        /// Assigns a 2-minute RT bin key for a given observed retention time.
        /// Matches the binning used by SSRCalc in PepAnalysisEngine.
        /// </summary>
        public static int GetRtBinKey(double observedRT)
            => (int)(2 * Math.Round(observedRT / 2d, 0));
    }
}
