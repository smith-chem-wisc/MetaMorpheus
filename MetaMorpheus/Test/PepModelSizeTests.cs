using EngineLayer;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// The PEP model is sized to its training set, and the PEP training q-value starts from a default per search type.
    /// </summary>
    [TestFixture]
    public static class PepModelSizeTests
    {
        [Test]
        public static void SmallTrainingSetsGetTheSmallModel()
        {
            var options = PepAnalysisEngine.GetTreeOptions(PepAnalysisEngine.SmallTrainingSetThreshold - 1, 42);
            Assert.That(options.NumberOfTrees, Is.EqualTo(50));
            Assert.That(options.NumberOfLeaves, Is.EqualTo(8));
            Assert.That(options.MinimumExampleCountPerLeaf, Is.EqualTo(20));
            Assert.That(options.LearningRate, Is.EqualTo(0.1));
        }

        [Test]
        public static void LargeTrainingSetsKeepTheEstablishedModel()
        {
            var options = PepAnalysisEngine.GetTreeOptions(PepAnalysisEngine.SmallTrainingSetThreshold, 42);
            Assert.That(options.NumberOfTrees, Is.EqualTo(400));
            Assert.That(options.NumberOfLeaves, Is.EqualTo(20));
            Assert.That(options.MinimumExampleCountPerLeaf, Is.EqualTo(10));
            Assert.That(options.LearningRate, Is.EqualTo(0.2));
        }

        [Test]
        [TestCase(false, false, 0.005)]
        [TestCase(false, true, 0.025)]
        [TestCase(true, false, 0.05)]
        [TestCase(true, true, 0.05)]
        public static void PepTrainingQValueDefaultsBySearchType(bool isRna, bool isTopDown, double expected)
        {
            Assert.That(PepAnalysisEngine.DefaultQValueCutoffForPepCalculation(isRna, isTopDown), Is.EqualTo(expected));
        }
    }
}
