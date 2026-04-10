using EngineLayer.FdrAnalysis;
using NUnit.Framework;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    [Category("RequiresKoinaAPI")]
    public class PrositRTIntegrationTests
    {
        [Test]
        [Explicit("Requires live Koina API access")]
        public void Prosit2019iRT_ReturnsValidPredictions_ForStandardPeptides()
        {
            // Standard peptide sequences valid for Prosit (length 7-30, standard AAs)
            var sequences = new List<string>
            {
                "PEPTIDER",
                "ELVISLIVESK",
                "AAADLETSSLK",
                "IGDYAGIK",
                "TASEFDSAIAQDK"
            };

            var model = new Prosit2019iRT();
            var inputs = sequences.Select(s => new RetentionTimePredictionInput(s)).ToList();
            var predictions = model.Predict(inputs);

            Assert.That(predictions, Has.Count.EqualTo(sequences.Count));
            foreach (var prediction in predictions)
            {
                Assert.That(prediction.PredictedRetentionTime, Is.Not.Null,
                    $"Prediction for {prediction.FullSequence} was null");
                Assert.That(prediction.PredictedRetentionTime.Value, Is.Not.NaN,
                    $"Prediction for {prediction.FullSequence} was NaN");
                Assert.That(prediction.PredictedRetentionTime.Value, Is.InRange(-200.0, 300.0),
                    $"Prediction for {prediction.FullSequence} was outside reasonable iRT range");
            }
        }

        [Test]
        [Explicit("Requires live Koina API access")]
        public void Prosit2019iRT_LinearCalibration_HasHighCorrelation()
        {
            // Generate 200 peptide-like sequences
            var rng = new Random(42);
            var aminoAcids = "ACDEFGHIKLMNPQRSTVWY";
            var sequences = new List<string>();
            for (int i = 0; i < 200; i++)
            {
                int len = rng.Next(7, 20);
                var seq = new string(Enumerable.Range(0, len)
                    .Select(_ => aminoAcids[rng.Next(aminoAcids.Length)])
                    .ToArray());
                sequences.Add(seq);
            }
            sequences = sequences.Distinct().ToList();

            var model = new Prosit2019iRT();
            var inputs = sequences.Select(s => new RetentionTimePredictionInput(s)).ToList();
            var predictions = model.Predict(inputs);

            // Filter to successful predictions
            var validPredictions = predictions
                .Where(p => p.PredictedRetentionTime.HasValue)
                .ToList();

            Assert.That(validPredictions.Count, Is.GreaterThanOrEqualTo(100),
                "Need at least 100 valid predictions for calibration test");

            // Simulate observed RTs as linear transform of iRT: observedRT = 1.5 * iRT + 10 + noise
            var rng2 = new Random(99);
            var iRTs = validPredictions.Select(p => p.PredictedRetentionTime.Value).ToList();
            var observedRTs = iRTs.Select(irt => 1.5 * irt + 10.0 + (rng2.NextDouble() - 0.5) * 2.0).ToList();

            (double slope, double intercept) = IrtCalibrationHelper.LinearRegression(iRTs, observedRTs);

            Assert.That(double.IsNaN(slope), Is.False, "Regression returned NaN slope");
            Assert.That(slope, Is.EqualTo(1.5).Within(1.5 * 0.10)); // within 10%
            Assert.That(intercept, Is.EqualTo(10.0).Within(10.0 * 0.20)); // within 20%
        }

        [Test]
        [Explicit("Requires live Koina API access")]
        public void Prosit2019iRT_WrappedInAdapter_ImplementsInterface()
        {
            var model = new Prosit2019iRT();
            IRetentionTimePredictor predictor = new RetentionTimeModelAdapter(model);

            Assert.That(predictor.IsIndexedRetentionTimeModel, Is.True);
            Assert.That(predictor.ModelName, Is.Not.Null.And.Not.Empty);

            var inputs = new List<RetentionTimePredictionInput>
            {
                new RetentionTimePredictionInput("PEPTIDER"),
                new RetentionTimePredictionInput("ELVISLIVESK")
            };

            var predictions = predictor.Predict(inputs);
            Assert.That(predictions, Has.Count.EqualTo(2));
            Assert.That(predictions[0].PredictedRetentionTime, Is.Not.Null);
        }
    }
}
