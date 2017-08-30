using EngineLayer.Gptmd;
using NUnit.Framework;
using Proteomics;

namespace Test
{
    [TestFixture]
    public static class ModFitsTest
    {
        #region Public Methods

        [Test]
        public static void TestModFits()
        {
            ModificationMotif.TryGetMotif("X", out ModificationMotif motif);
            ModificationWithMass attemptToLocalize = new ModificationWithMass(null, null, motif, TerminusLocalization.Any, double.NaN);

            Protein protein = new Protein("M", null);
            int peptideOneBasedIndex = 1;
            int peptideLength = 1;
            int proteinOneBasedIndex = 1;
            Assert.IsTrue(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("M", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, TerminusLocalization.Any, double.NaN);
            Assert.IsTrue(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("N", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, TerminusLocalization.Any, double.NaN);
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("Mx", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, TerminusLocalization.Any, double.NaN);
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("Mr", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, TerminusLocalization.Any, double.NaN);
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("xM", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, TerminusLocalization.Any, double.NaN);
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("Nxs", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, TerminusLocalization.Any, double.NaN);

            protein = new Protein("MNRS", null);
            peptideOneBasedIndex = 1;
            peptideLength = 1;
            proteinOneBasedIndex = 1;
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("Nxs", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, TerminusLocalization.Any, double.NaN);

            protein = new Protein("MNRS", null);
            peptideOneBasedIndex = 1;
            peptideLength = 1;
            proteinOneBasedIndex = 1;
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));
            peptideOneBasedIndex = 2;
            peptideLength = 1;
            proteinOneBasedIndex = 2;
            Assert.IsTrue(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));
            protein = new Protein("MNRN", null);
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));
        }

        #endregion Public Methods
    }
}