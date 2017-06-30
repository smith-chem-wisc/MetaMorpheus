using EngineLayer;
using MzLibUtil;
using NUnit.Framework;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class SearchModesTest
    {

        #region Public Methods

        [Test]
        public static void TestSearchModeTest()
        {
            MassDiffAcceptor sm = new TestSearchMode("My custom");
            Assert.IsTrue(sm.Accepts(2, 2) >= 0);
            Assert.IsTrue(sm.Accepts(0.5, 4) >= 0);
            Assert.IsFalse(sm.Accepts(0.5, 0.5) >= 0);
            Assert.IsTrue(sm.Accepts(1, 1) >= 0);
            Assert.AreEqual(2, sm.GetAllowedPrecursorMassIntervals(0.5).First().allowedInterval.Minimum);
            Assert.AreEqual(0.5, sm.GetAllowedPrecursorMassIntervals(2).First().allowedInterval.Minimum);
        }

        [Test]
        public static void TestDotSearchMode()
        {
            var dsm1 = new DotMassDiffAcceptor("test1", new double[] { 0, 1 }, new Tolerance(ToleranceUnit.Absolute, 0.1));

            Assert.IsTrue(dsm1.Accepts(1000, 1000) >= 0);
            Assert.IsTrue(dsm1.Accepts(1000, 1000 + 0.1 / 2) >= 0);
            Assert.IsFalse(dsm1.Accepts(1000, 1000 + 0.1 * 2) >= 0);
            Assert.IsTrue(dsm1.Accepts(1000 + 0.1 / 2, 1000) >= 0);
            Assert.IsFalse(dsm1.Accepts(1000 + 0.1 * 2, 1000) >= 0);

            Assert.IsTrue(dsm1.Accepts(1000 + 1, 1000 + 0.1 / 2) >= 0);
            Assert.IsFalse(dsm1.Accepts(1000 + 1, 1000 + 0.1 * 2) >= 0);
            Assert.IsTrue(dsm1.Accepts(1000 + 1 + 0.1 / 2, 1000) >= 0);
            Assert.IsFalse(dsm1.Accepts(1000 + 1 + 0.1 * 2, 1000) >= 0);

            var theList = dsm1.GetAllowedPrecursorMassIntervals(100).ToList();

            Assert.AreEqual(99.9, theList[0].allowedInterval.Minimum);
            Assert.AreEqual(100.1, theList[0].allowedInterval.Maximum);
            Assert.AreEqual(100.9, theList[1].allowedInterval.Minimum);
            Assert.AreEqual(101.1, theList[1].allowedInterval.Maximum);

            var dsm2 = new DotMassDiffAcceptor("test2", new double[] { 0, 1 }, new Tolerance(ToleranceUnit.PPM, 5));

            Assert.IsTrue(dsm2.Accepts(1000, 1000) >= 0);

            Assert.IsTrue(dsm2.Accepts(1000 * (1 + 5.0 / 1e6 / 1.0000001), 1000) >= 0); // FIRST VARIES WITHIN 5 PPM OF SECOND
            Assert.IsTrue(dsm2.Accepts(1000 * (1 - 5.0 / 1e6 / 1.0000001), 1000) >= 0); // FIRST VARIES WITHIN 5 PPM OF SECOND

            Assert.IsFalse(dsm2.Accepts(1000, 1000 * (1 - 5.0 / 1e6 / 1.0000001)) >= 0); // VERY CAREFUL

            Assert.IsFalse(dsm2.Accepts(1000 * (1 + 5.0 / 1e6 * 1.0000001), 1000) >= 0); // FIRST VARIES WITHIN 5 PPM OF SECOND
            Assert.IsFalse(dsm2.Accepts(1000 * (1 - 5.0 / 1e6 * 1.0000001), 1000) >= 0); // FIRST VARIES WITHIN 5 PPM OF SECOND

            Assert.IsTrue(dsm2.Accepts(1000, 1000 * (1 + 5.0 / 1e6 * 1.0000001)) >= 0); // VERY CAREFUL

            var theList2 = dsm2.GetAllowedPrecursorMassIntervals(1000).ToList();

            Assert.IsTrue(theList2[0].allowedInterval.Contains(1000));

            Assert.IsTrue(1000 * (1 + 5.0 / 1e6 / 1.0000001) < theList2[0].allowedInterval.Maximum);
            Assert.IsTrue(1000 * (1 - 5.0 / 1e6 / 1.0000001) > theList2[0].allowedInterval.Minimum);
            Assert.IsTrue(1000 * (1 + 5.0 / 1e6 * 1.0000001) > theList2[0].allowedInterval.Maximum);
            Assert.IsTrue(1000 * (1 - 5.0 / 1e6 * 1.0000001) < theList2[0].allowedInterval.Minimum);

            Assert.IsTrue(theList2[1].allowedInterval.Contains(1001));
        }

        [Test]
        public static void TestIntervalsSearchMode()
        {
        }

        #endregion Public Methods

        #region Private Classes

        // Accept if scanPrecursorMass*peptideMass>=1.
        private class TestSearchMode : MassDiffAcceptor
        {

            #region Public Constructors

            public TestSearchMode(string fileNameAddition) : base(fileNameAddition)
            {
            }

            #endregion Public Constructors

            #region Public Methods

            public override int Accepts(double scanPrecursorMass, double peptideMass)
            {
                return scanPrecursorMass * peptideMass >= 1 ? 1 : -1;
            }

            public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
            {
                yield return new AllowedIntervalWithNotch(new DoubleRange(1 / peptideMonoisotopicMass, double.MaxValue), 1);
            }

            #endregion Public Methods

        }

        #endregion Private Classes

    }
}