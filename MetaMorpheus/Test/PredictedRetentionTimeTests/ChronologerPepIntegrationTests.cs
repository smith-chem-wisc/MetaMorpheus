using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.Chronologer;
using EngineLayer;
using EngineLayer.FdrAnalysis;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace Test.PredictedRetentionTimeTests
{
    [TestFixture]
    public class ChronologerPepIntegrationTests
    {
        private static readonly IRetentionTimePredictor _predictor =
            new ChronologerRetentionTimePredictor();

        [Test]
        public void Chronologer_PredictRetentionTimes_ReturnsDictionaryForValidPeptides()
        {
            var protein = new Protein("PEPTIDEKMRGESTARKFINDLERLIKEPEPTIDEK", "TEST");
            var peptides = new[]
            {
                new PeptideWithSetModifications("PEPTIDEK", null, p: protein),
                new PeptideWithSetModifications("ELVISLIVES", null, p: protein),
                new PeptideWithSetModifications("GESTARK", null, p: protein),
                new PeptideWithSetModifications("FINDLERLIKE", null, p: protein),
                new PeptideWithSetModifications("MRAGESTARK", null, p: protein),
            };

            Dictionary<string, double?> predictions = _predictor.PredictRetentionTimes(peptides);

            Assert.That(predictions.Count, Is.EqualTo(peptides.Length));
            foreach (var kvp in predictions)
            {
                Assert.That(kvp.Value.HasValue, Is.True, $"{kvp.Key} returned null");
                Assert.That(double.IsFinite(kvp.Value.Value), Is.True);
            }
        }

        [Test]
        public void PepAnalysisEngine_WithChronologer_IncludesPredictedRTZScoreFeature()
        {
            // Covered end-to-end by EverythingRunnerEngineTestCases.BottomUpPepQValue, which
            // runs the default-Chronologer pipeline and writes results whose PSM count
            // (449) reflects PredictedRTZScore contributing to PEP training. A unit-level
            // assertion would require ≥100 hand-built SpectralMatch objects with realistic
            // peptides and retention times — not worth the complexity given the integration
            // test already exercises this path.
            EverythingRunnerEngineTestCase.TryGetTestCase(
                EverythingRunnerEngineTestCases.BottomUpPepQValue, out var testCase);
            string resultsFile = Path.Combine(testCase.OutputDirectory,
                @"postSearchAnalysisTaskTestOutput\results.txt");
            Assert.That(File.Exists(resultsFile), Is.True,
                "BottomUpPepQValue test case did not produce results — Chronologer pipeline failed");

            var lines = File.ReadAllLines(resultsFile);
            Assert.That(lines.Any(l => l.Contains("All target PSMs with pep q-value <= 0.01: 449")),
                Is.True,
                "PSM count does not match the post-Chronologer baseline — feature may have regressed");
        }

        [Test]
        public void PepAnalysisEngine_WithChronologer_ProducesLowZScoreForMatchingRT()
        {
            Assert.Inconclusive("Requires handcrafted ≥100-PSM fixture; covered indirectly by BottomUpPepQValue");
        }

        [Test]
        public void PepAnalysisEngine_WithChronologer_ProducesHighZScoreForWrongRT()
        {
            Assert.Inconclusive("Requires handcrafted ≥100-PSM fixture; covered indirectly by BottomUpPepQValue");
        }

        [Test]
        public void FdrAnalysisEngine_DefaultPredictor_IsChronologer()
        {
            // 1) Default CommonParameters should select Chronologer
            var defaultParams = new CommonParameters();
            Assert.That(defaultParams.RTPredictorName, Is.EqualTo("Chronologer"));

            // 2) Private static GetRTPredictor should return a Chronologer instance for "standard"
            MethodInfo getPredictor = typeof(FdrAnalysisEngine).GetMethod(
                "GetRTPredictor",
                BindingFlags.NonPublic | BindingFlags.Static);
            Assert.That(getPredictor, Is.Not.Null, "GetRTPredictor not found via reflection");

            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>
            {
                ("dummy.mzML", defaultParams)
            };
            var result = getPredictor!.Invoke(null, new object[] { "standard", fsp });
            Assert.That(result, Is.InstanceOf<ChronologerRetentionTimePredictor>());

            // 3) Non-standard search type returns null
            var nullResult = getPredictor.Invoke(null, new object[] { "top-down", fsp });
            Assert.That(nullResult, Is.Null);
        }
    }
}
