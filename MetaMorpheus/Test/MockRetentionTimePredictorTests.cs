using EngineLayer.FdrAnalysis;
using NUnit.Framework;
using PredictionClients.Koina.AbstractClasses;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Http;

namespace Test
{
    /// <summary>
    /// Returns deterministic iRT predictions from a pre-loaded dictionary.
    /// Sequences not in the dictionary get null predictions (as if the model rejected them).
    /// </summary>
    public class MockRetentionTimePredictor : IRetentionTimePredictor
    {
        private readonly Dictionary<string, double> _predictions;
        public bool IsIndexedRetentionTimeModel => true;
        public string ModelName => "MockRTModel";
        public bool WasCalled { get; private set; }
        public int CallCount { get; private set; }

        public MockRetentionTimePredictor(Dictionary<string, double> predictions)
            => _predictions = predictions;

        public List<PeptideRTPrediction> Predict(List<RetentionTimePredictionInput> inputs)
        {
            WasCalled = true;
            CallCount++;
            return inputs.Select(i => new PeptideRTPrediction(
                FullSequence: i.FullSequence,
                ValidatedFullSequence: i.FullSequence,
                PredictedRetentionTime: _predictions.TryGetValue(i.FullSequence, out var irt) ? irt : null,
                IsIndexed: true,
                Warning: null
            )).ToList();
        }
    }

    /// <summary>
    /// Always throws HttpRequestException to simulate a network failure.
    /// </summary>
    public class FailingRetentionTimePredictor : IRetentionTimePredictor
    {
        public bool IsIndexedRetentionTimeModel => true;
        public string ModelName => "FailingRTModel";

        public List<PeptideRTPrediction> Predict(List<RetentionTimePredictionInput> inputs)
            => throw new HttpRequestException("Simulated network failure");
    }

    [TestFixture]
    public class MockRetentionTimePredictorTests
    {
        // Test 1: PredictedRTZScore IS in standard trainingInfos
        [Test]
        public void PsmData_StandardTrainingInfos_ContainsPredictedRTZScore()
        {
            Assert.That(PsmData.trainingInfos["standard"], Does.Contain("PredictedRTZScore"));
        }

        // Test 2: PredictedRTZScore is NOT in top-down, crosslink, RNA trainingInfos
        [Test]
        public void PsmData_NonStandardTrainingInfos_DoNotContainPredictedRTZScore()
        {
            Assert.That(PsmData.trainingInfos["top-down"], Does.Not.Contain("PredictedRTZScore"));
            Assert.That(PsmData.trainingInfos["crosslink"], Does.Not.Contain("PredictedRTZScore"));
            Assert.That(PsmData.trainingInfos["RNA"], Does.Not.Contain("PredictedRTZScore"));
        }

        // Test 3: assumedAttributeDirection for PredictedRTZScore is -1
        // (lower z-score = better RT match = more likely true positive)
        [Test]
        public void PsmData_AssumedAttributeDirection_PredictedRTZScore_IsNegativeOne()
        {
            Assert.That(PsmData.assumedAttributeDirection["PredictedRTZScore"], Is.EqualTo(-1));
        }

        // Test 4: PepAnalysisEngine calls the mock predictor exactly once
        // (even if many PSMs share the same sequence — deduplication)
        [Test]
        public void PepAnalysisEngine_WithMockPredictor_CallsPredictorExactlyOnce()
        {
            // TODO: Requires PepAnalysisEngine constructor to accept IRetentionTimePredictor (Step 7)
            // and ComputePredictedRTValues implementation (Step 8)
            Assert.Inconclusive("Requires Steps 7-8 implementation");
        }

        // Test 5: PepAnalysisEngine with mock predictor includes PredictedRTZScore in TrainingVariables
        [Test]
        public void PepAnalysisEngine_WithMockPredictor_IncludesFeatureInTrainingVariables()
        {
            // TODO: Requires PepAnalysisEngine constructor to accept IRetentionTimePredictor (Step 7)
            Assert.Inconclusive("Requires Steps 7-8 implementation");
        }

        // Test 6: PepAnalysisEngine with FailingRetentionTimePredictor does NOT throw
        [Test]
        public void PepAnalysisEngine_WithFailingPredictor_DoesNotThrow()
        {
            // TODO: Requires fallback wiring (Step 9)
            Assert.Inconclusive("Requires Step 9 implementation");
        }

        // Test 7: PepAnalysisEngine with FailingRetentionTimePredictor removes feature from TrainingVariables
        [Test]
        public void PepAnalysisEngine_WithFailingPredictor_RemovesFeatureFromTrainingVariables()
        {
            // TODO: Requires fallback wiring (Step 9)
            Assert.Inconclusive("Requires Step 9 implementation");
        }

        // Test 8: PepAnalysisEngine with null predictor never includes feature
        [Test]
        public void PepAnalysisEngine_WithNullPredictor_DoesNotIncludeFeature()
        {
            // TODO: Requires PepAnalysisEngine constructor changes (Step 7)
            Assert.Inconclusive("Requires Step 7 implementation");
        }

        // Test 9: When <100 unique calibration peptides available, feature is skipped for that file
        [Test]
        public void PepAnalysisEngine_InsufficientCalibrationPeptides_SkipsFeatureForFile()
        {
            // TODO: Requires ComputePredictedRTValues implementation (Step 8)
            Assert.Inconclusive("Requires Step 8 implementation");
        }

        // Test 10: PsmData.PredictedRTZScore is low (< 1.0) when prediction matches observation
        [Test]
        public void PsmData_PredictedRTZScore_IsLow_WhenPredictionMatchesObservation()
        {
            // TODO: Requires GetPredictedRTZScore implementation (Step 10)
            Assert.Inconclusive("Requires Step 10 implementation");
        }

        // Test 11: PsmData.PredictedRTZScore is maxed (10.0) when prediction is wildly wrong
        [Test]
        public void PsmData_PredictedRTZScore_IsMaxed_WhenPredictionIsWrong()
        {
            // TODO: Requires GetPredictedRTZScore implementation (Step 10)
            Assert.Inconclusive("Requires Step 10 implementation");
        }

        // Test 12: PsmData.PredictedRTZScore is sentinel (10.0) when sequence has no prediction
        [Test]
        public void PsmData_PredictedRTZScore_IsSentinel_WhenSequenceNotPredicted()
        {
            // TODO: Requires GetPredictedRTZScore implementation (Step 10)
            Assert.Inconclusive("Requires Step 10 implementation");
        }
    }
}
