using NUnit.Framework;
using EngineLayer.FdrAnalysis;
using Chromatography.RetentionTimePrediction;
using Chromatography;
using System.Collections.Generic;

namespace Test.PredictedRetentionTimeTests
{
    // ── Mock predictors ───────────────────────────────────────────────────

    /// <summary>
    /// Returns a fixed prediction for every sequence. Used to test the happy path.
    /// </summary>
    internal class ConstantRetentionTimePredictor : IRetentionTimePredictor
    {
        private readonly double _value;
        public string PredictorName => "Constant";
        public SeparationType SeparationType => SeparationType.HPLC;

        public ConstantRetentionTimePredictor(double value) => _value = value;

        public double? PredictRetentionTime(IRetentionPredictable peptide,
            out RetentionTimeFailureReason? failureReason)
        {
            failureReason = null;
            return _value;
        }

        public string? GetFormattedSequence(IRetentionPredictable peptide,
            out RetentionTimeFailureReason? failureReason)
        {
            failureReason = null;
            return peptide.BaseSequence;
        }
    }

    /// <summary>
    /// Always throws. Used to test fallback behavior.
    /// </summary>
    internal class FailingRetentionTimePredictor : IRetentionTimePredictor
    {
        public string PredictorName => "Failing";
        public SeparationType SeparationType => SeparationType.HPLC;

        public double? PredictRetentionTime(IRetentionPredictable peptide,
            out RetentionTimeFailureReason? failureReason)
            => throw new System.Exception("Simulated predictor failure");

        public Dictionary<string, double?> PredictRetentionTimes(
            IEnumerable<IRetentionPredictable> peptides)
            => throw new System.Exception("Simulated batch predictor failure");

        public string? GetFormattedSequence(IRetentionPredictable peptide,
            out RetentionTimeFailureReason? failureReason)
        {
            failureReason = null;
            return null;
        }
    }

    [TestFixture]
    public class MockPredictorTests
    {
        // Test 1: PredictedRTZScore is in standard trainingInfos
        [Test]
        public void PsmData_StandardTrainingInfos_ContainsPredictedRTZScore()
        {
            Assert.That(PsmData.trainingInfos["standard"],
                Does.Contain("PredictedRTZScore"));
        }

        // Test 2: PredictedRTZScore has direction -1
        [Test]
        public void PsmData_PredictedRTZScore_HasNegativeDirection()
        {
            Assert.That(PsmData.assumedAttributeDirection["PredictedRTZScore"],
                Is.EqualTo(-1));
        }

        // Test 3: PredictedRTZScore is NOT in non-standard trainingInfos
        [Test]
        public void PsmData_NonStandardTrainingInfos_DoNotContainPredictedRTZScore()
        {
            foreach (var key in PsmData.trainingInfos.Keys)
            {
                if (key == "standard") continue;
                Assert.That(PsmData.trainingInfos[key],
                    Does.Not.Contain("PredictedRTZScore"),
                    $"Search type '{key}' should not contain PredictedRTZScore");
            }
        }

        // Test 4: null predictor -> feature not in TrainingVariables
        [Test]
        public void PepAnalysisEngine_NullPredictor_DoesNotIncludeFeature()
        {
            // Construct a minimal PepAnalysisEngine with null predictor.
            // Ask developer how to construct with minimal PSMs for this test.
            // TrainingVariables should not contain "PredictedRTZScore".
            Assert.Inconclusive("Requires SpectralMatch construction — implement after Step 6");
        }

        // Test 5: FailingRetentionTimePredictor -> no crash, feature removed
        [Test]
        public void PepAnalysisEngine_FailingPredictor_DoesNotThrow()
        {
            Assert.Inconclusive("Requires SpectralMatch construction — implement after Step 6");
        }

        // Test 6: FailingRetentionTimePredictor -> PredictedRTZScore removed from TrainingVariables
        [Test]
        public void PepAnalysisEngine_FailingPredictor_RemovesFeatureFromTrainingVariables()
        {
            Assert.Inconclusive("Requires SpectralMatch construction — implement after Step 6");
        }

        // Test 7: <100 calibration peptides -> sentinel for that file
        [Test]
        public void PepAnalysisEngine_InsufficientCalibrationPeptides_UsesSentinel()
        {
            Assert.Inconclusive("Requires SpectralMatch construction — implement after Step 6");
        }
    }
}
