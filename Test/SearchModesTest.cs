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
            Assert.IsTrue(sm.Accepts(2, 2) >= 0);
            Assert.IsTrue(sm.Accepts(0.5, 4) >= 0);
            Assert.IsFalse(sm.Accepts(0.5, 0.5) >= 0);
            Assert.IsTrue(sm.Accepts(1, 1) >= 0);
            Assert.AreEqual(2, sm.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(0.5).First().AllowedInterval.Minimum);
            Assert.AreEqual(0.5, sm.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(2).First().AllowedInterval.Minimum);
        }

        [Test]
        public static void TestDotSearchMode()
        {
            var dsm1 = new DotMassDiffAcceptor("test1", new double[] { 0, 1 }, new AbsoluteTolerance(0.1));

            Assert.IsTrue(dsm1.Accepts(1000, 1000) >= 0);
            Assert.IsTrue(dsm1.Accepts(1000, 1000 + 0.1 / 2) >= 0);
            Assert.IsFalse(dsm1.Accepts(1000, 1000 + 0.1 * 2) >= 0);
            Assert.IsTrue(dsm1.Accepts(1000 + 0.1 / 2, 1000) >= 0);
            Assert.IsFalse(dsm1.Accepts(1000 + 0.1 * 2, 1000) >= 0);

            Assert.IsTrue(dsm1.Accepts(1000 + 1, 1000 + 0.1 / 2) >= 0);
            Assert.IsFalse(dsm1.Accepts(1000 + 1, 1000 + 0.1 * 2) >= 0);
            Assert.IsTrue(dsm1.Accepts(1000 + 1 + 0.1 / 2, 1000) >= 0);
            Assert.IsFalse(dsm1.Accepts(1000 + 1 + 0.1 * 2, 1000) >= 0);

            var theList = dsm1.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(100).ToList();

            Assert.AreEqual(99.9, theList[0].AllowedInterval.Minimum);
            Assert.AreEqual(100.1, theList[0].AllowedInterval.Maximum);
            Assert.AreEqual(100.9, theList[1].AllowedInterval.Minimum);
            Assert.AreEqual(101.1, theList[1].AllowedInterval.Maximum);

            var dsm2 = new DotMassDiffAcceptor("test2", new double[] { 0, 1 }, new PpmTolerance(5));

            Assert.IsTrue(dsm2.Accepts(1000, 1000) >= 0);

            Assert.IsTrue(dsm2.Accepts(1000 * (1 + 5.0 / 1e6 / 1.0000001), 1000) >= 0); // FIRST VARIES WITHIN 5 PPM OF SECOND
            Assert.IsTrue(dsm2.Accepts(1000 * (1 - 5.0 / 1e6 / 1.0000001), 1000) >= 0); // FIRST VARIES WITHIN 5 PPM OF SECOND

            Assert.IsFalse(dsm2.Accepts(1000, 1000 * (1 - 5.0 / 1e6 / 1.0000001)) >= 0); // VERY CAREFUL

            Assert.IsFalse(dsm2.Accepts(1000 * (1 + 5.0 / 1e6 * 1.0000001), 1000) >= 0); // FIRST VARIES WITHIN 5 PPM OF SECOND
            Assert.IsFalse(dsm2.Accepts(1000 * (1 - 5.0 / 1e6 * 1.0000001), 1000) >= 0); // FIRST VARIES WITHIN 5 PPM OF SECOND

            Assert.IsTrue(dsm2.Accepts(1000, 1000 * (1 + 5.0 / 1e6 * 1.0000001)) >= 0); // VERY CAREFUL

            var theList2 = dsm2.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(1000).ToList();

            Assert.IsTrue(theList2[0].AllowedInterval.Contains(1000));

            Assert.IsTrue(1000 * (1 + 5.0 / 1e6 / 1.0000001) < theList2[0].AllowedInterval.Maximum);
            Assert.IsTrue(1000 * (1 - 5.0 / 1e6 / 1.0000001) > theList2[0].AllowedInterval.Minimum);
            Assert.IsTrue(1000 * (1 + 5.0 / 1e6 * 1.0000001) > theList2[0].AllowedInterval.Maximum);
            Assert.IsTrue(1000 * (1 - 5.0 / 1e6 * 1.0000001) < theList2[0].AllowedInterval.Minimum);

            Assert.IsTrue(theList2[1].AllowedInterval.Contains(1001));
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