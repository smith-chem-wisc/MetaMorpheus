using FragmentGeneration;
using NUnit.Framework;

namespace Test
{
    [TestFixture]
    public class Class1
    {
        [Test]
        public void TestExclusions()
        {
            double tolExclude = 0.0075;
            double[] exclude = new double[5] { 1, 2, 3, 4, 5 };

            Assert.IsTrue(Exclusions.DoNotExclude(0, tolExclude, exclude));
            Assert.IsFalse(Exclusions.DoNotExclude(1, tolExclude, exclude));

            exclude = Exclusions.PopulateExcludeList();

            // A real reverse phosphorylation sequence
            Assert.IsFalse(Exclusions.DoNotExclude(-76.134779, tolExclude, exclude));

            // Do not exclude zero
            Assert.IsTrue(Exclusions.DoNotExclude(0, tolExclude, exclude));

            // Glycine - do not exclude
            Assert.IsTrue(Exclusions.DoNotExclude(57.02146, tolExclude, exclude));

            // Valine - exclude
            Assert.IsFalse(Exclusions.DoNotExclude(99.06841, tolExclude, exclude));

            // Combo - exclude
            Assert.IsFalse(Exclusions.DoNotExclude(57.02146 + 99.06841, tolExclude, exclude));

            // Lysine - do not exclude
            Assert.IsTrue(Exclusions.DoNotExclude(128.09496, tolExclude, exclude));

            // Lysine combo - do not exclude
            Assert.IsTrue(Exclusions.DoNotExclude(128.09496 + 57.02146, tolExclude, exclude));

            // Lysine combo - do not exclude
            Assert.IsTrue(Exclusions.DoNotExclude(128.09496 + 99.06841, tolExclude, exclude));
        }
    }
}