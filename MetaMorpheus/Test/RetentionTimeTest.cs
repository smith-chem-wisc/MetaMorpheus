using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Proteomics.RetentionTimePrediction;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    public sealed class RetentionTimeTest
    {
        [Test]
        public static void SimpleReversePhaseRetentionTimeTest()
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
                double actual = calc.ScoreSequence(peptide);

                Assert.That(expected, Is.EqualTo(actual).Within(0.01));
            }
        }
    }
}