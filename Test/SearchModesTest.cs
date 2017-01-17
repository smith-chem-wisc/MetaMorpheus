using InternalLogicEngineLayer;
using NUnit.Framework;
using Spectra;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class SearchModesTest
    {

        #region Public Methods

        [Test]
        public void TestSearchModeTest()
        {
            SearchMode sm = new TestSearchMode("My custom");
            Assert.IsTrue(sm.Accepts(2, 2));
            Assert.IsTrue(sm.Accepts(0.5, 4));
            Assert.IsFalse(sm.Accepts(0.5, 0.5));
            Assert.IsTrue(sm.Accepts(1, 1));
            Assert.AreEqual(2, sm.GetAllowedPrecursorMassIntervals(0.5).First().Minimum);
            Assert.AreEqual(0.5, sm.GetAllowedPrecursorMassIntervals(2).First().Minimum);
        }

        [Test]
        public void TestDotSearchMode()
        {
            var dsm1 = new DotSearchMode("test1", new double[] { 0, 1 }, new Tolerance(ToleranceUnit.Absolute, 0.1));

            Assert.IsTrue(dsm1.Accepts(1000, 1000));
            Assert.IsTrue(dsm1.Accepts(1000, 1000 + 0.1 / 2));
            Assert.IsFalse(dsm1.Accepts(1000, 1000 + 0.1 * 2));
            Assert.IsTrue(dsm1.Accepts(1000 + 0.1 / 2, 1000));
            Assert.IsFalse(dsm1.Accepts(1000 + 0.1 * 2, 1000));

            Assert.IsTrue(dsm1.Accepts(1000 + 1, 1000 + 0.1 / 2));
            Assert.IsFalse(dsm1.Accepts(1000 + 1, 1000 + 0.1 * 2));
            Assert.IsTrue(dsm1.Accepts(1000 + 1 + 0.1 / 2, 1000));
            Assert.IsFalse(dsm1.Accepts(1000 + 1 + 0.1 * 2, 1000));

            var theList = dsm1.GetAllowedPrecursorMassIntervals(100).ToList();

            Assert.AreEqual(99.9, theList[0].Minimum);
            Assert.AreEqual(100.1, theList[0].Maximum);
            Assert.AreEqual(100.9, theList[1].Minimum);
            Assert.AreEqual(101.1, theList[1].Maximum);

            var dsm2 = new DotSearchMode("test2", new double[] { 0, 1 }, new Tolerance(ToleranceUnit.PPM, 5));

            Assert.IsTrue(dsm2.Accepts(1000, 1000));

            Assert.IsTrue(dsm2.Accepts(1000 * (1 + 5.0 / 1e6 / 1.0000001), 1000)); // FIRST VARIES WITHIN 5 PPM OF SECOND
            Assert.IsTrue(dsm2.Accepts(1000 * (1 - 5.0 / 1e6 / 1.0000001), 1000)); // FIRST VARIES WITHIN 5 PPM OF SECOND

            Assert.IsFalse(dsm2.Accepts(1000, 1000 * (1 - 5.0 / 1e6 / 1.0000001))); // VERY CAREFUL

            Assert.IsFalse(dsm2.Accepts(1000 * (1 + 5.0 / 1e6 * 1.0000001), 1000)); // FIRST VARIES WITHIN 5 PPM OF SECOND
            Assert.IsFalse(dsm2.Accepts(1000 * (1 - 5.0 / 1e6 * 1.0000001), 1000)); // FIRST VARIES WITHIN 5 PPM OF SECOND

            Assert.IsTrue(dsm2.Accepts(1000, 1000 * (1 + 5.0 / 1e6 * 1.0000001))); // VERY CAREFUL

            var theList2 = dsm2.GetAllowedPrecursorMassIntervals(1000).ToList();

            Assert.IsTrue(theList2[0].Contains(1000));

            Assert.IsTrue(1000 * (1 + 5.0 / 1e6 / 1.0000001) < theList2[0].Maximum);
            Assert.IsTrue(1000 * (1 - 5.0 / 1e6 / 1.0000001) > theList2[0].Minimum);
            Assert.IsTrue(1000 * (1 + 5.0 / 1e6 * 1.0000001) > theList2[0].Maximum);
            Assert.IsTrue(1000 * (1 - 5.0 / 1e6 * 1.0000001) < theList2[0].Minimum);

            Assert.IsTrue(theList2[1].Contains(1001));
        }

        [Test]
        public void TestIntervalsSearchMode()
        {
        }

        #endregion Public Methods

        #region Private Classes

        // Accept if scanPrecursorMass*peptideMass>=1.
        private class TestSearchMode : SearchMode
        {

            #region Public Constructors

            public TestSearchMode(string fileNameAddition) : base(fileNameAddition)
            {
            }

            #endregion Public Constructors

            #region Public Methods

            public override bool Accepts(double scanPrecursorMass, double peptideMass)
            {
                return scanPrecursorMass * peptideMass >= 1;
            }

            public override IEnumerable<DoubleRange> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
            {
                yield return new DoubleRange(1 / peptideMonoisotopicMass, double.MaxValue);
            }

            public override string SearchModeString()
            {
                return "scanPrecursorMass * peptideMass >= 1";
            }

            #endregion Public Methods

        }

        #endregion Private Classes

    }
}