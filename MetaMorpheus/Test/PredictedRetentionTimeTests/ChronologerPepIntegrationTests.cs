using NUnit.Framework;
using EngineLayer.FdrAnalysis;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.Chronologer;

namespace Test.PredictedRetentionTimeTests
{
    [TestFixture]
    public class ChronologerPepIntegrationTests
    {
        [Test]
        public void Chronologer_PredictRetentionTimes_ReturnsDictionaryForValidPeptides()
        {
            Assert.Inconclusive("Implement in Step 8");
        }

        [Test]
        public void PepAnalysisEngine_WithChronologer_IncludesPredictedRTZScoreFeature()
        {
            Assert.Inconclusive("Implement in Step 8");
        }

        [Test]
        public void PepAnalysisEngine_WithChronologer_ProducesLowZScoreForMatchingRT()
        {
            Assert.Inconclusive("Implement in Step 8");
        }

        [Test]
        public void PepAnalysisEngine_WithChronologer_ProducesHighZScoreForWrongRT()
        {
            Assert.Inconclusive("Implement in Step 8");
        }

        [Test]
        public void FdrAnalysisEngine_DefaultPredictor_IsChronologer()
        {
            Assert.Inconclusive("Implement in Step 8");
        }
    }
}
