using EngineLayer;
using NUnit.Framework;
using System;
using System.Threading.Tasks;

namespace Test
{
    [TestFixture]
    public static class TailorScoreTest
    {
        /// <summary>
        /// Builds a distribution of <paramref name="lowCount"/> observations at
        /// <paramref name="low"/> and <paramref name="highCount"/> at <paramref name="high"/>.
        /// </summary>
        private static ScoreDistribution Distribution(int lowCount, double low, int highCount, double high)
        {
            var d = new ScoreDistribution();
            for (int i = 0; i < lowCount; i++) d.Add(low);
            for (int i = 0; i < highCount; i++) d.Add(high);
            return d;
        }

        [Test]
        public static void PercentileIsRecoveredToWithinHalfABin()
        {
            // 900 candidates at 1.0 and 100 at 10.0: the 990th of 1000 sorted ascending is 10.0.
            var d = Distribution(900, 1.0, 100, 10.0);
            Assert.That(d.Count, Is.EqualTo(1000));
            Assert.That(d.TryGetPercentile(0.99, out double q99), Is.True);
            Assert.That(q99, Is.EqualTo(10.0).Within(d.BinWidth / 2));

            // 990 at 1.0 and 10 at 20.0 puts the same rank down at 1.0 instead.
            var d2 = Distribution(990, 1.0, 10, 20.0);
            Assert.That(d2.TryGetPercentile(0.99, out double q99Low), Is.True);
            Assert.That(q99Low, Is.EqualTo(1.0).Within(d2.BinWidth / 2));
        }

        [Test]
        public static void PercentileTracksTheRequestedRank()
        {
            var d = Distribution(500, 2.0, 500, 8.0);
            Assert.That(d.TryGetPercentile(0.25, out double q25), Is.True);
            Assert.That(q25, Is.EqualTo(2.0).Within(d.BinWidth / 2));
            Assert.That(d.TryGetPercentile(0.75, out double q75), Is.True);
            Assert.That(q75, Is.EqualTo(8.0).Within(d.BinWidth / 2));
        }

        [Test]
        public static void EmptyDistributionHasNoPercentile()
        {
            var d = new ScoreDistribution();
            Assert.That(d.Count, Is.EqualTo(0));
            Assert.That(d.TryGetPercentile(0.99, out _), Is.False);
        }

        [Test]
        public static void ScoresAboveTheClampAreRefusedRatherThanUnderstated()
        {
            // Every observation lands in the clamped top bin, where the true value is known only
            // to be at or above the threshold. A number here would silently understate the divisor.
            var d = new ScoreDistribution();
            double aboveClamp = d.ClampThreshold + 100;
            for (int i = 0; i < 1000; i++) d.Add(aboveClamp);

            Assert.That(d.Count, Is.EqualTo(1000));
            Assert.That(d.TryGetPercentile(0.99, out _), Is.False);
            Assert.That(TailorScoreCalculator.Calculate(aboveClamp, d), Is.NaN);
        }

        [Test]
        public static void NonFiniteAndNegativeScoresAreIgnored()
        {
            var d = new ScoreDistribution();
            d.Add(double.NaN);
            d.Add(-1.0);
            Assert.That(d.Count, Is.EqualTo(0));

            d.Add(5.0);
            Assert.That(d.Count, Is.EqualTo(1));
        }

        [Test]
        public static void ConcurrentAddsAreNotLost()
        {
            // The search loop is protein-major and parallel, so one scan's histogram is written
            // from every thread; a lost update would bias the percentile.
            var d = new ScoreDistribution();
            const int perThread = 10_000;
            const int threads = 8;

            Parallel.For(0, threads, _ =>
            {
                for (int i = 0; i < perThread; i++) d.Add(3.0);
            });

            Assert.That(d.Count, Is.EqualTo(perThread * threads));
            Assert.That(d.TryGetPercentile(0.99, out double q99), Is.True);
            Assert.That(q99, Is.EqualTo(3.0).Within(d.BinWidth / 2));
        }

        [Test]
        public static void CalibratedScoreDividesByTheHighPercentileNotTheLowOne()
        {
            // The regression this test exists for: dividing by the *bottom* of the distribution,
            // floored at 1, returns the raw score unchanged for any realistic candidate set.
            var d = Distribution(900, 1.0, 100, 10.0);
            double rawScore = 30.0;

            double tailor = TailorScoreCalculator.Calculate(rawScore, d);

            Assert.That(tailor, Is.Not.NaN);
            Assert.That(tailor, Is.EqualTo(rawScore / 10.125).Within(0.01));
            Assert.That(tailor, Is.LessThan(rawScore / 2),
                "dividing by a high percentile must move the score substantially");
        }

        [Test]
        public static void SpectraWithTooFewCandidatesAreNotCalibrated()
        {
            var justUnder = Distribution(TailorScoreCalculator.MinimumCandidates - 1, 1.0, 0, 0);
            Assert.That(TailorScoreCalculator.Calculate(30.0, justUnder), Is.NaN);

            var atThreshold = Distribution(TailorScoreCalculator.MinimumCandidates, 1.0, 0, 0);
            Assert.That(TailorScoreCalculator.Calculate(30.0, atThreshold), Is.Not.NaN);
        }

        [Test]
        public static void NullDistributionIsNotCalibrated()
        {
            Assert.That(TailorScoreCalculator.Calculate(30.0, null), Is.NaN);
        }

        [Test]
        public static void CalibrationIsMonotonicInTheRawScore()
        {
            var d = Distribution(900, 1.0, 100, 10.0);
            double low = TailorScoreCalculator.Calculate(10.0, d);
            double high = TailorScoreCalculator.Calculate(20.0, d);
            Assert.That(high, Is.GreaterThan(low));
        }

        [Test]
        public static void BadConstructorAndPercentileArgumentsThrow()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new ScoreDistribution(binWidth: 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => new ScoreDistribution(binCount: 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => new ScoreDistribution().TryGetPercentile(0, out _));
            Assert.Throws<ArgumentOutOfRangeException>(() => new ScoreDistribution().TryGetPercentile(1.5, out _));
        }
    }
}
