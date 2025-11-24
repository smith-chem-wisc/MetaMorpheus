using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using Chromatography.RetentionTimePrediction.SSRCalc;
using Omics.Modifications;
using Chromatography.RetentionTimePrediction.Chronologer;

namespace Test
{
    [TestFixture]
    public sealed class RetentionTimeTest
    {
        [Test]
        public static void SimpleReversePhaseRetentionTimeTest_SSRCalc_Direct()
        {
            Dictionary<string, double> peptidesAndHyrophobicities = new Dictionary<string, double>
            {
                { "QSHFANAEPEQK", 11.27 },
                { "SDLFENLQNYR", 30.44 },
                { "SLPSEFEPINLR", 33.12 }
            };

            SSRCalc3 calc = new SSRCalc3("SSRCalc 3.0 (300A)", SSRCalc3.Column.A300);

            foreach (string peptideSequence in peptidesAndHyrophobicities.Keys)
            {
                var peptide = new PeptideWithSetModifications(peptideSequence, new Dictionary<string, Modification>());
                double expected = peptidesAndHyrophobicities[peptideSequence];
                double actual = calc.ScoreSequence(peptide.BaseSequence);

                Assert.That(expected, Is.EqualTo(actual).Within(0.01));
            }
        }

        [Test]
        public static void SimpleReversePhaseRetentionTimeTest_SSRCalc_Structure()
        {
            Dictionary<string, double> peptidesAndHyrophobicities = new Dictionary<string, double>
            {
                { "QSHFANAEPEQK", 11.27 },
                { "SDLFENLQNYR", 30.44 },
                { "SLPSEFEPINLR", 33.12 }
            };

            var predictor = new SSRCalc3RetentionTimePredictor();

            foreach (string peptideSequence in peptidesAndHyrophobicities.Keys)
            {
                var peptide = new PeptideWithSetModifications(peptideSequence, new Dictionary<string, Modification>());
                double expected = peptidesAndHyrophobicities[peptideSequence];
                double? actual = predictor.PredictRetentionTime(peptide, out var failureReason);

                Assert.That(failureReason, Is.Null);
                Assert.That(expected, Is.EqualTo(actual).Within(0.01));
            }
        }

        [Test]
        public static void SimpleReversePhaseRetentionTimeTest_Chronologer_Structure()
        {
            Dictionary<string, double> peptidesAndHyrophobicities = new Dictionary<string, double>
            {
                { "QSHFANAEPEQK", 2.75 },
                { "SDLFENLQNYR", 14.81 },
                { "SLPSEFEPINLR", 15.88 }
            };

            var predictor = new ChronologerRetentionTimePredictor();

            foreach (string peptideSequence in peptidesAndHyrophobicities.Keys)
            {
                var peptide = new PeptideWithSetModifications(peptideSequence, new Dictionary<string, Modification>());
                double expected = peptidesAndHyrophobicities[peptideSequence];
                double? actual = predictor.PredictRetentionTime(peptide, out var failureReason);

                Assert.That(failureReason, Is.Null);
                Assert.That(expected, Is.EqualTo(actual).Within(0.01));
            }
        }
    }
}