using EngineLayer;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public static class SearchModesTest
    {
        [Test]
        public static void TestSearchModeTest()
        {
            MassDiffAcceptor sm = new TestSearchMode("My custom");
            Assert.That(sm.Accepts(2, 2) >= 0);
            Assert.That(sm.Accepts(0.5, 4) >= 0);
            Assert.That(!(sm.Accepts(0.5, 0.5) >= 0));
            Assert.That(sm.Accepts(1, 1) >= 0);
            Assert.That(sm.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(0.5).First().AllowedInterval.Minimum, Is.EqualTo(2));
            Assert.That(sm.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(2).First().AllowedInterval.Minimum, Is.EqualTo(0.5));
        }

        [Test]
        public static void TestDotSearchMode()
        {
            var dsm1 = new DotMassDiffAcceptor("test1", new double[] { 0, 1 }, new AbsoluteTolerance(0.1));

            Assert.That(dsm1.Accepts(1000, 1000) >= 0);
            Assert.That(dsm1.Accepts(1000, 1000 + 0.1 / 2) >= 0);
            Assert.That(!(dsm1.Accepts(1000, 1000 + 0.1 * 2) >= 0));
            Assert.That(dsm1.Accepts(1000 + 0.1 / 2, 1000) >= 0);
            Assert.That(!(dsm1.Accepts(1000 + 0.1 * 2, 1000) >= 0));

            Assert.That(dsm1.Accepts(1000 + 1, 1000 + 0.1 / 2) >= 0);
            Assert.That(!(dsm1.Accepts(1000 + 1, 1000 + 0.1 * 2) >= 0));
            Assert.That(dsm1.Accepts(1000 + 1 + 0.1 / 2, 1000) >= 0);
            Assert.That(!(dsm1.Accepts(1000 + 1 + 0.1 * 2, 1000) >= 0));

            var theList = dsm1.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(100).ToList();

            Assert.That(theList[0].AllowedInterval.Minimum, Is.EqualTo(99.9));
            Assert.That(theList[0].AllowedInterval.Maximum, Is.EqualTo(100.1));
            Assert.That(theList[1].AllowedInterval.Minimum, Is.EqualTo(100.9));
            Assert.That(theList[1].AllowedInterval.Maximum, Is.EqualTo(101.1));

            var dsm2 = new DotMassDiffAcceptor("test2", new double[] { 0, 1 }, new PpmTolerance(5));

            Assert.That(dsm2.Accepts(1000, 1000) >= 0);

            Assert.That(dsm2.Accepts(1000 * (1 + 5.0 / 1e6 / 1.0000001), 1000) >= 0); // FIRST VARIES WITHIN 5 PPM OF SECOND
            Assert.That(dsm2.Accepts(1000 * (1 - 5.0 / 1e6 / 1.0000001), 1000) >= 0); // FIRST VARIES WITHIN 5 PPM OF SECOND

            Assert.That(!(dsm2.Accepts(1000, 1000 * (1 - 5.0 / 1e6 / 1.0000001)) >= 0)); // VERY CAREFUL

            Assert.That(!(dsm2.Accepts(1000 * (1 + 5.0 / 1e6 * 1.0000001), 1000) >= 0)); // FIRST VARIES WITHIN 5 PPM OF SECOND
            Assert.That(!(dsm2.Accepts(1000 * (1 - 5.0 / 1e6 * 1.0000001), 1000) >= 0)); // FIRST VARIES WITHIN 5 PPM OF SECOND

            Assert.That(dsm2.Accepts(1000, 1000 * (1 + 5.0 / 1e6 * 1.0000001)) >= 0); // VERY CAREFUL

            var theList2 = dsm2.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(1000).ToList();

            Assert.That(theList2[0].AllowedInterval.Contains(1000));

            Assert.That(1000 * (1 + 5.0 / 1e6 / 1.0000001) < theList2[0].AllowedInterval.Maximum);
            Assert.That(1000 * (1 - 5.0 / 1e6 / 1.0000001) > theList2[0].AllowedInterval.Minimum);
            Assert.That(1000 * (1 + 5.0 / 1e6 * 1.0000001) > theList2[0].AllowedInterval.Maximum);
            Assert.That(1000 * (1 - 5.0 / 1e6 * 1.0000001) < theList2[0].AllowedInterval.Minimum);

            Assert.That(theList2[1].AllowedInterval.Contains(1001));
        }

        // Accept if scanPrecursorMass*peptideMass>=1.
        private class TestSearchMode : MassDiffAcceptor
        {
            public TestSearchMode(string fileNameAddition) : base(fileNameAddition)
            {
            }

            public override int Accepts(double scanPrecursorMass, double peptideMass)
            {
                return scanPrecursorMass * peptideMass >= 1 ? 1 : -1;
            }

            public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
            {
                yield return new AllowedIntervalWithNotch(new DoubleRange(1 / peptideMonoisotopicMass, double.MaxValue), 1);
            }

            public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
            {
                yield return new AllowedIntervalWithNotch(new DoubleRange(double.MinValue, 1 / peptideMonoisotopicMass), 1);
            }

            public override string ToProseString()
            {
                throw new NotImplementedException();
            }
        }
    }
}