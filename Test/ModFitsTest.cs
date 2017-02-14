using EngineLayer.Gptmd;
using NUnit.Framework;
using Proteomics;

namespace Test
{
    [TestFixture]
    public class ModFitsTest
    {

        #region Public Methods

        [Test]
        public static void TestModFits()
        {
            ModificationMotif motif;
            ModificationMotif.TryGetMotif("X", out motif);
            ModificationWithMass attemptToLocalize = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, double.NaN, null, null, null);
            string proteinBaseSequence = "M";
            int peptideOneBasedIndex = 1;
            int peptideLength = 1;
            int proteinOneBasedIndex = 1;
            int proteinLength = 1;
            Assert.IsTrue(GptmdEngine.ModFits(attemptToLocalize, proteinBaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex, proteinLength));

            ModificationMotif.TryGetMotif("M", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, double.NaN, null, null, null);
            Assert.IsTrue(GptmdEngine.ModFits(attemptToLocalize, proteinBaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex, proteinLength));

            ModificationMotif.TryGetMotif("N", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, double.NaN, null, null, null);
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, proteinBaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex, proteinLength));

            ModificationMotif.TryGetMotif("Mx", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, double.NaN, null, null, null);
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, proteinBaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex, proteinLength));

            ModificationMotif.TryGetMotif("Mr", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, double.NaN, null, null, null);
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, proteinBaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex, proteinLength));

            ModificationMotif.TryGetMotif("xM", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, double.NaN, null, null, null);
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, proteinBaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex, proteinLength));

            ModificationMotif.TryGetMotif("Nxs", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, double.NaN, null, null, null);
            proteinBaseSequence = "MNRS";
            peptideOneBasedIndex = 1;
            peptideLength = 1;
            proteinOneBasedIndex = 1;
            proteinLength = 4;
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, proteinBaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex, proteinLength));

            ModificationMotif.TryGetMotif("Nxs", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, double.NaN, null, null, null);
            proteinBaseSequence = "MNRS";
            peptideOneBasedIndex = 1;
            peptideLength = 1;
            proteinOneBasedIndex = 1;
            proteinLength = 4;
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, proteinBaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex, proteinLength));
            peptideOneBasedIndex = 2;
            peptideLength = 1;
            proteinOneBasedIndex = 2;
            proteinLength = 4;
            Assert.IsTrue(GptmdEngine.ModFits(attemptToLocalize, proteinBaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex, proteinLength));
            proteinBaseSequence = "MNRN";
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, proteinBaseSequence, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex, proteinLength));
        }

        #endregion Public Methods

    }
}