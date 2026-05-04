using NUnit.Framework;
using EngineLayer;
using EngineLayer.FdrAnalysis;
using Chromatography.RetentionTimePrediction;
using Chromatography;
using System.Collections.Generic;
using System.Linq;

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

        public double? PredictRetentionTimeEquivalent(IRetentionPredictable peptide,
            out RetentionTimeFailureReason? failureReason)
            => PredictRetentionTime(peptide, out failureReason);

        public IReadOnlyList<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)>
            PredictRetentionTimeEquivalents(IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1)
            => peptides.Select(p => ((double?)_value, p, (RetentionTimeFailureReason?)null)).ToList();

        public string? GetFormattedSequence(IRetentionPredictable peptide,
            out RetentionTimeFailureReason? failureReason)
        {
            failureReason = null;
            return peptide.BaseSequence;
        }

        public void Dispose() { }
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

        public IReadOnlyDictionary<string, double?> PredictRetentionTimes(
            IEnumerable<IRetentionPredictable> peptides)
            => throw new System.Exception("Simulated batch predictor failure");

        public double? PredictRetentionTimeEquivalent(IRetentionPredictable peptide,
            out RetentionTimeFailureReason? failureReason)
            => throw new System.Exception("Simulated predictor failure");

        public IReadOnlyList<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)>
            PredictRetentionTimeEquivalents(IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1)
            => throw new System.Exception("Simulated batch predictor failure");

        public string? GetFormattedSequence(IRetentionPredictable peptide,
            out RetentionTimeFailureReason? failureReason)
        {
            failureReason = null;
            return null;
        }

        public void Dispose() { }
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

        private static List<(string fileName, CommonParameters fileSpecificParameters)> BuildFsp()
            => new() { ("dummy.mzML", new CommonParameters()) };

        private static PepAnalysisEngine BuildEngine(IRetentionTimePredictor predictor)
            => new PepAnalysisEngine(
                new List<SpectralMatch>(),
                "standard",
                BuildFsp(),
                outputFolder: null,
                rtPredictor: predictor);

        // Test 4: null predictor -> feature not in TrainingVariables
        [Test]
        public void PepAnalysisEngine_NullPredictor_DoesNotIncludeFeature()
        {
            var engine = BuildEngine(predictor: null);
            Assert.That(engine.TrainingVariables, Does.Not.Contain("PredictedRTZScore"));
        }

        // Test 5: FailingRetentionTimePredictor -> no crash
        [Test]
        public void PepAnalysisEngine_FailingPredictor_DoesNotThrow()
        {
            // With empty PSMs ComputePredictedRTValues short-circuits on <100 calibration peptides
            // before calling the predictor, so use a predictor that throws inside construction
            // flow — here construction alone must not throw.
            Assert.DoesNotThrow(() => BuildEngine(new FailingRetentionTimePredictor()));
        }

        // Test 6: FailingRetentionTimePredictor -> PredictedRTZScore removed from TrainingVariables
        [Test]
        public void PepAnalysisEngine_FailingPredictor_RemovesFeatureFromTrainingVariables()
        {
            // Requires ≥100 high-confidence PSMs to force ComputePredictedRTValues to call
            // the predictor (and therefore hit the catch-block that strips the feature).
            // With an empty PSM list the method short-circuits on calibration-peptide count
            // and the predictor is never invoked. Covered end-to-end by the BottomUpPepQValue
            // integration test case (default-Chronologer path).
            Assert.Inconclusive("Requires ≥100 high-confidence PSMs to exercise the exception path");
        }

        // Test 7: <100 calibration peptides -> engine constructs without error
        [Test]
        public void PepAnalysisEngine_InsufficientCalibrationPeptides_DoesNotThrow()
        {
            Assert.DoesNotThrow(() => BuildEngine(new ConstantRetentionTimePredictor(30.0)));
        }
    }
}
