using System;
using NUnit.Framework;
using EngineLayer.FdrAnalysis;
using System.Collections.Generic;
using System.Linq;

namespace Test.PredictedRetentionTimeTests
{
    [TestFixture]
    public class IrtCalibrationTests
    {
        // ── OLS Regression ────────────────────────────────────────────────

        [Test]
        public void OlsRegression_PerfectLinearData_ReturnsExactSlopeAndIntercept()
        {
            // y = 2x + 5
            var x = new List<double> { 1, 2, 3, 4, 5 };
            var y = x.Select(v => 2 * v + 5).ToList();

            var (slope, intercept) = IrtCalibrationHelper.OlsRegression(x, y);

            Assert.That(slope, Is.EqualTo(2.0).Within(1e-9));
            Assert.That(intercept, Is.EqualTo(5.0).Within(1e-9));
        }

        [Test]
        public void OlsRegression_NoisyData_ReturnsReasonableFit()
        {
            // y ≈ 1.5x + 10 with small noise
            var x = Enumerable.Range(1, 50).Select(i => (double)i).ToList();
            var rng = new System.Random(42);
            var y = x.Select(v => 1.5 * v + 10 + (rng.NextDouble() - 0.5)).ToList();

            var (slope, intercept) = IrtCalibrationHelper.OlsRegression(x, y);

            Assert.That(slope, Is.EqualTo(1.5).Within(0.1));
            Assert.That(intercept, Is.EqualTo(10.0).Within(1.0));
        }

        [Test]
        public void OlsRegression_InsufficientData_ReturnsNaN()
        {
            var x = new List<double> { 1.0 };
            var y = new List<double> { 2.0 };

            var (slope, intercept) = IrtCalibrationHelper.OlsRegression(x, y);

            Assert.That(double.IsNaN(slope), Is.True);
            Assert.That(double.IsNaN(intercept), Is.True);
        }

        [Test]
        public void OlsRegression_EmptyData_ReturnsNaN()
        {
            var (slope, intercept) = IrtCalibrationHelper.OlsRegression(
                new List<double>(), new List<double>());

            Assert.That(double.IsNaN(slope), Is.True);
            Assert.That(double.IsNaN(intercept), Is.True);
        }

        // ── Bin Assignment ────────────────────────────────────────────────

        [Test]
        public void GetTwoMinuteBin_CorrectlyBinsRetentionTimes()
        {
            // Bin width is 2 minutes; bin = (int)(rt / 2)
            Assert.That(IrtCalibrationHelper.GetTwoMinuteBin(0.0), Is.EqualTo(0));
            Assert.That(IrtCalibrationHelper.GetTwoMinuteBin(1.9), Is.EqualTo(0));
            Assert.That(IrtCalibrationHelper.GetTwoMinuteBin(2.0), Is.EqualTo(1));
            Assert.That(IrtCalibrationHelper.GetTwoMinuteBin(3.9), Is.EqualTo(1));
            Assert.That(IrtCalibrationHelper.GetTwoMinuteBin(30.0), Is.EqualTo(15));
        }

        // ── Bin Statistics ────────────────────────────────────────────────

        [Test]
        public void ComputeBinStatistics_SingleBin_ReturnsCorrectMeanAndStdDev()
        {
            // All residuals in the same 2-minute bin, mean=0, stddev=1
            var residuals = new List<(double observedRT, double residual)>
            {
                (5.0, -1.0),
                (5.1,  0.0),
                (5.2,  1.0),
                (5.3,  0.0),
                (5.4,  0.0),
            };

            var binStats = IrtCalibrationHelper.ComputeBinStatistics(residuals);

            int bin = IrtCalibrationHelper.GetTwoMinuteBin(5.0);
            Assert.That(binStats.ContainsKey(bin), Is.True);
            Assert.That(binStats[bin].Item1, Is.EqualTo(0.0).Within(1e-9)); // mean
            Assert.That(binStats[bin].Item2, Is.GreaterThan(0)); // stddev
        }

        [Test]
        public void ComputeBinStatistics_EmptyInput_ReturnsEmptyDictionary()
        {
            var binStats = IrtCalibrationHelper.ComputeBinStatistics(
                new List<(double, double)>());

            Assert.That(binStats, Is.Empty);
        }

        // ── StdDev Stabilization ──────────────────────────────────────────

        [Test]
        public void StabilizeBinStDevs_TooTightBins_AreRaisedToMinimum()
        {
            var binStats = new Dictionary<int, Tuple<double, double>>
            {
                { 0, Tuple.Create(0.0, 0.001) }, // very tight — below minimum
                { 1, Tuple.Create(0.0, 2.5) },   // normal
            };

            IrtCalibrationHelper.StabilizeBinStDevs(binStats);

            // Tight bin should be raised; normal bin should be unchanged
            Assert.That(binStats[0].Item2, Is.GreaterThan(0.001));
            Assert.That(binStats[1].Item2, Is.EqualTo(2.5).Within(1e-9));
        }

        // ── Z-Score Computation ───────────────────────────────────────────

        [Test]
        public void ComputeZScore_ResidualWithinOneSigma_ReturnsLessThanOne()
        {
            double residual = 0.5;
            double mean = 0.0;
            double stddev = 1.0;

            double zScore = IrtCalibrationHelper.ComputeZScore(residual, mean, stddev);

            Assert.That(zScore, Is.LessThan(1.0));
        }

        [Test]
        public void ComputeZScore_ResidualFarFromMean_ReturnsHighValue()
        {
            double residual = 100.0;
            double mean = 0.0;
            double stddev = 1.0;

            double zScore = IrtCalibrationHelper.ComputeZScore(residual, mean, stddev);

            Assert.That(zScore, Is.GreaterThan(5.0));
        }

        [Test]
        public void ComputeZScore_ZeroStdDev_DoesNotThrow()
        {
            Assert.DoesNotThrow(() =>
                IrtCalibrationHelper.ComputeZScore(1.0, 0.0, 0.0));
        }

        [Test]
        public void ComputeZScore_IsCappedAtSentinelValue()
        {
            // Z-scores should be capped at 10.0 (the sentinel)
            double zScore = IrtCalibrationHelper.ComputeZScore(1000.0, 0.0, 0.001);

            Assert.That(zScore, Is.EqualTo(10.0));
        }
    }
}
