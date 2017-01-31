using EngineLayer;
using EngineLayer.Gptmd;
using MassSpectrometry;
using NUnit.Framework;

using Spectra;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class GptmdEngineTest
    {

        #region Public Methods

        [Test]
        public static void TestGptmdEngine()
        {
            List<NewPsmWithFdr> allResultingIdentifications = null;
            var gptmdModifications = new List<MetaMorpheusModification> { new MetaMorpheusModification("name", ModificationType.AminoAcidResidue, 'N', null, '\0', double.NaN, double.NaN, 21.981943, new Chemistry.ChemicalFormula("H-1 Na1")) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>>();
            Tolerance precursorMassTolerance = new Tolerance(ToleranceUnit.PPM, 10);
            bool isotopeErrors = false;

            allResultingIdentifications = new List<NewPsmWithFdr>();
            var engine = new GptmdEngine(allResultingIdentifications, isotopeErrors, gptmdModifications, combos, precursorMassTolerance);
            var res = (GptmdResults)engine.Run();
            Assert.AreEqual(0, res.Mods.Count);

            PsmParent newPsm = new TestParentSpectrumMatch(588.22520189093 + 21.981943);
            var parentProtein = new Protein("NNNNN", "accession", new Dictionary<int, List<MetaMorpheusModification>>(), null, null, null, null, null, 0, false, false);
            var modPep = new PeptideWithPossibleModifications(1, 5, parentProtein, 0, "ugh");
            //var twoBasedVariableAndLocalizeableModificationss = new Dictionary<int, MorpheusModification>();
            List<MetaMorpheusModification> variableModifications = new List<MetaMorpheusModification>();
            var peptidesWithSetModifications = new HashSet<PeptideWithSetModifications> { modPep.GetPeptideWithSetModifications(variableModifications, 4096, 3, new List<MetaMorpheusModification>()).First() };
            Tolerance fragmentTolerance = null;
            IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile = null;
            var thisPSM = new PsmWithMultiplePossiblePeptides(newPsm, peptidesWithSetModifications, fragmentTolerance, myMsDataFile, new List<ProductType> { ProductType.B, ProductType.Y });
            var newPsmWithFDR = new NewPsmWithFdr(thisPSM, 1, 0, 0);
            allResultingIdentifications.Add(newPsmWithFDR);

            engine = new GptmdEngine(allResultingIdentifications, isotopeErrors, gptmdModifications, combos, precursorMassTolerance);
            res = (GptmdResults)engine.Run();
            Assert.AreEqual(1, res.Mods.Count);
            Assert.AreEqual(5, res.Mods["accession"].Count);
        }

        #endregion Public Methods

        #region Private Classes

        private class TestParentSpectrumMatch : PsmParent
        {

            #region Public Constructors

            public TestParentSpectrumMatch(double scanPrecursorMass) : base(null, double.NaN, double.NaN, scanPrecursorMass, 0, 0, 0, double.NaN, double.NaN, double.NaN)
            {
            }

            #endregion Public Constructors

            #region Public Methods

            public override CompactPeptide GetCompactPeptide(List<MetaMorpheusModification> variableModifications, List<MetaMorpheusModification> localizeableModifications, List<MetaMorpheusModification> fixedModifications)
            {
                throw new NotImplementedException();
            }

            #endregion Public Methods

        }

        #endregion Private Classes

    }
}