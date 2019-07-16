using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Proteomics.RetentionTimePrediction;
using System;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    public sealed class RetentionTimeTest
    {
        private static object[,] _peptides300A;

        [Test]
        public static void SimpleReversePhaseRetentionTimeTest()
        {
            _peptides300A = new object[,]
                {
                {"QSHFANAEPEQK", 11.27},
                {"SDLFENLQNYR", 30.44},
                {"SLPSEFEPINLR", 33.12}
                };

            SSRCalc3 calc = new SSRCalc3("SSRCalc 3.0 (300A)", SSRCalc3.Column.A300);

            for (int i = 0; i < _peptides300A.GetLength(0); i++)
            {
                var peptide = new PeptideWithSetModifications((string)_peptides300A[i, 0], new Dictionary<string, Modification>());
                double expected = (double)_peptides300A[i, 1];
                double actual = calc.ScoreSequence(peptide);

                // Round the returned value to match the values presented
                // in the supporting information of the SSRCalc 3 publication.
                // First cast to float, since the added precision of the double
                // caused the double representation of 12.805 to round to 12.80
                // instead of 12.81.  When diffed with 12.81 the result was
                // 0.005000000000002558.
                double actualRound = Math.Round((float)actual, 2);

                // Extra conditional added to improve debugging of issues.
                if (Math.Abs(expected - actual) > 0.005)
                {
                    Assert.AreEqual(expected, actualRound, "Peptide {0}", peptide);
                }
            }
        }
    }
}