using EngineLayer.FdrAnalysis;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class IrtCalibrationTests
    {
        // Test 1: Perfect linear relationship
        // Given: points exactly on y = 2x + 5
        // Assert: slope ≈ 2.0, intercept ≈ 5.0 (within 1e-6)
        [Test]
        public void LinearRegression_PerfectLine_ReturnsExactSlopeAndIntercept()
        {
            // Generate 200 points on y = 2x + 5
            var x = Enumerable.Range(0, 200).Select(i => (double)i).ToList();
            var y = x.Select(xi => 2.0 * xi + 5.0).ToList();

            (double slope, double intercept) = IrtCalibrationHelper.LinearRegression(x, y);

            Assert.That(slope, Is.EqualTo(2.0).Within(1e-6));
            Assert.That(intercept, Is.EqualTo(5.0).Within(1e-6));
        }

        // Test 2: Near-perfect real-world-like data
        // Given: 200 synthetic (iRT, observedRT) pairs where observedRT = 1.5*iRT + 10 + small_noise
        // Assert: slope within 5% of 1.5, intercept within 5% of 10, R² > 0.99
        [Test]
        public void LinearRegression_RealWorldData_HasHighCorrelation()
        {
            var rng = new Random(42);
            var x = new List<double>();
            var y = new List<double>();
            for (int i = 0; i < 200; i++)
            {
                double xi = i * 0.5; // iRT values from 0 to 99.5
                double noise = (rng.NextDouble() - 0.5) * 0.5; // small noise ±0.25
                x.Add(xi);
                y.Add(1.5 * xi + 10.0 + noise);
            }

            (double slope, double intercept) = IrtCalibrationHelper.LinearRegression(x, y);

            Assert.That(slope, Is.EqualTo(1.5).Within(1.5 * 0.05)); // within 5%
            Assert.That(intercept, Is.EqualTo(10.0).Within(10.0 * 0.05)); // within 5%

            // Verify R² > 0.99 by computing it from residuals
            double yMean = y.Average();
            double ssTot = y.Sum(yi => (yi - yMean) * (yi - yMean));
            double ssRes = x.Zip(y, (xi, yi) => yi - (slope * xi + intercept))
                            .Sum(r => r * r);
            double rSquared = 1.0 - ssRes / ssTot;
            Assert.That(rSquared, Is.GreaterThan(0.99));
        }

        // Test 3: Insufficient data returns NaN sentinel
        // Given: only 5 data points (below the 100 minimum)
        [Test]
        public void LinearRegression_TooFewPoints_ReturnsNaNSentinel()
        {
            var x = new List<double> { 1, 2, 3, 4, 5 };
            var y = new List<double> { 2, 4, 6, 8, 10 };

            (double slope, double intercept) = IrtCalibrationHelper.LinearRegression(x, y);

            Assert.That(double.IsNaN(slope), Is.True);
            Assert.That(double.IsNaN(intercept), Is.True);
        }

        // Test 4: Bin statistics on clean data
        // Given: residuals drawn from N(0,1) assigned to a single bin
        [Test]
        public void ComputeBinStatistics_NormalResiduals_HasCorrectMeanAndStdDev()
        {
            // Generate ~1000 residuals approximately from N(0,1) using Box-Muller
            var rng = new Random(123);
            var residuals = new List<double>();
            for (int i = 0; i < 1000; i++)
            {
                double u1 = rng.NextDouble();
                double u2 = rng.NextDouble();
                double z = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
                residuals.Add(z);
            }

            var residualsByBin = new Dictionary<int, List<double>>
            {
                { 10, residuals }
            };

            var stats = IrtCalibrationHelper.ComputeBinStatistics(residualsByBin);

            Assert.That(stats.ContainsKey(10), Is.True);
            Assert.That(stats[10].Item1, Is.EqualTo(0.0).Within(0.1)); // mean ≈ 0
            Assert.That(stats[10].Item2, Is.EqualTo(1.0).Within(0.2)); // stddev ≈ 1
        }

        // Test 5: StdDev stabilization replaces outlier bins (too tight)
        // Given: one bin with stddev = 0.001, globalStddev = 1.0
        [Test]
        public void StabilizeBinStDevs_ReplacesTooTightBin_WithGlobalStdDev()
        {
            var stats = new Dictionary<int, Tuple<double, double>>
            {
                { 10, Tuple.Create(0.0, 0.001) } // mean=0, stddev=0.001 (too tight)
            };

            IrtCalibrationHelper.StabilizeBinStDevs(stats, globalStdDev: 1.0);

            Assert.That(stats[10].Item2, Is.EqualTo(1.0));
        }

        // Test 6: StdDev stabilization replaces too-wide bins
        // Given: one bin with stddev = 5.0, globalStddev = 1.0 (ratio > 3)
        [Test]
        public void StabilizeBinStDevs_ReplacesTooWideBin_WithGlobalStdDev()
        {
            var stats = new Dictionary<int, Tuple<double, double>>
            {
                { 10, Tuple.Create(0.0, 5.0) } // mean=0, stddev=5.0 (too wide, ratio=5)
            };

            IrtCalibrationHelper.StabilizeBinStDevs(stats, globalStdDev: 1.0);

            Assert.That(stats[10].Item2, Is.EqualTo(1.0));
        }

        // Test 7: Z-score computation — perfect prediction
        // residual=0, binMean=0, binStddev=1 → zScore = 0.0
        [Test]
        public void ComputeZScore_PerfectPrediction_ReturnsZero()
        {
            double zScore = IrtCalibrationHelper.ComputeZScore(
                residual: 0.0, binMean: 0.0, binStdDev: 1.0);

            Assert.That(zScore, Is.EqualTo(0.0).Within(1e-10));
        }

        // Test 8: Z-score computation — large deviation is clamped
        // residual=999, binMean=0, binStddev=1 → zScore = 10.0 (maxZScore clamp)
        [Test]
        public void ComputeZScore_LargeDeviation_IsClamped()
        {
            double zScore = IrtCalibrationHelper.ComputeZScore(
                residual: 999.0, binMean: 0.0, binStdDev: 1.0);

            Assert.That(zScore, Is.EqualTo(IrtCalibrationHelper.MaxZScore));
        }

        // Test 9: Z-score computation — NaN stddev is handled
        // stddev=0 → zScore = 10.0 (not infinity, not NaN)
        [Test]
        public void ComputeZScore_ZeroStdDev_ReturnsSentinel()
        {
            double zScore = IrtCalibrationHelper.ComputeZScore(
                residual: 1.0, binMean: 0.0, binStdDev: 0.0);

            Assert.That(zScore, Is.EqualTo(IrtCalibrationHelper.MaxZScore));
        }

        // Test 10: RT binning uses 2-minute windows
        [Test]
        public void RtBinning_UsesCorrectTwoMinuteWindows()
        {
            // observedRT=5.3 → (int)(2 * Math.Round(5.3/2, 0)) = (int)(2 * Math.Round(2.65, 0)) = (int)(2*3) = 6
            Assert.That(IrtCalibrationHelper.GetRtBinKey(5.3), Is.EqualTo(6));

            // observedRT=4.7 → (int)(2 * Math.Round(4.7/2, 0)) = (int)(2 * Math.Round(2.35, 0)) = (int)(2*2) = 4
            Assert.That(IrtCalibrationHelper.GetRtBinKey(4.7), Is.EqualTo(4));

            // observedRT=0.0 → bin 0
            Assert.That(IrtCalibrationHelper.GetRtBinKey(0.0), Is.EqualTo(0));

            // observedRT=1.0 → (int)(2 * Math.Round(0.5, 0)) = (int)(2*0) = 0
            // (banker's rounding: 0.5 rounds to 0)
            Assert.That(IrtCalibrationHelper.GetRtBinKey(1.0), Is.EqualTo(0));

            // observedRT=60.0 → (int)(2 * Math.Round(30.0, 0)) = 60
            Assert.That(IrtCalibrationHelper.GetRtBinKey(60.0), Is.EqualTo(60));
        }
    }
}
