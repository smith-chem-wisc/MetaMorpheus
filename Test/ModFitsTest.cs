using EngineLayer.Gptmd;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;

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
            ModificationWithMass attemptToLocalize = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, new List<double> { double.NaN }, null, null, null);

            Protein protein = new Protein("M", null, null, null, null, null, null, null, null, false, false, null);
            int peptideOneBasedIndex = 1;
            int peptideLength = 1;
            int proteinOneBasedIndex = 1;
            Assert.IsTrue(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("M", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, new List<double> { double.NaN }, null, null, null);
            Assert.IsTrue(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("N", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, new List<double> { double.NaN }, null, null, null);
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("Mx", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, new List<double> { double.NaN }, null, null, null);
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("Mr", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, new List<double> { double.NaN }, null, null, null);
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("xM", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, new List<double> { double.NaN }, null, null, null);
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("Nxs", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, new List<double> { double.NaN }, null, null, null);

            protein = new Protein("MNRS", null, null, null, null, null, null, null, null, false, false, null);
            peptideOneBasedIndex = 1;
            peptideLength = 1;
            proteinOneBasedIndex = 1;
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));

            ModificationMotif.TryGetMotif("Nxs", out motif);
            attemptToLocalize = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, new List<double> { double.NaN }, null, null, null);

            protein = new Protein("MNRS", null, null, null, null, null, null, null, null, false, false, null);
            peptideOneBasedIndex = 1;
            peptideLength = 1;
            proteinOneBasedIndex = 1;
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));
            peptideOneBasedIndex = 2;
            peptideLength = 1;
            proteinOneBasedIndex = 2;
            Assert.IsTrue(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));
            protein = new Protein("MNRN", null, null, null, null, null, null, null, null, false, false, null);
            Assert.IsFalse(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));
        }

        #endregion Public Methods

    }
}