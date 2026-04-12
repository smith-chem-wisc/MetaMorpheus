using NUnit.Framework;

namespace Test.PredictedRetentionTimeTests
{
    [TestFixture]
    [Category("RequiresKoinaAPI")]
    public class KoinaPepIntegrationTests
    {
        [Test]
        [Explicit("Requires live Koina API access")]
        public void PepAnalysisEngine_WithProsit2019iRT_IncludesPredictedRTZScoreFeature()
        {
            Assert.Inconclusive("Implement after Koina PR merges or alongside");
        }

        [Test]
        [Explicit("Requires live Koina API access")]
        public void Prosit2019iRT_PredictRetentionTimes_ReturnsValidIRTValues()
        {
            Assert.Inconclusive("Implement after Koina PR merges or alongside");
        }
    }
}
