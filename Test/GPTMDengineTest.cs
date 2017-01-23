using InternalLogicEngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using OldInternalLogic;
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
            var gptmdModifications = new List<MorpheusModification> { new MorpheusModification("name", ModificationType.AminoAcidResidue, 'N', 42, null, null, '\0', double.NaN, false, null) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>>();
            double tol = 0.1;
            bool isotopeErrors = false;

            allResultingIdentifications = new List<NewPsmWithFdr>();
            var engine = new GptmdEngine(allResultingIdentifications, isotopeErrors, gptmdModifications, combos, tol);
            var res = (GptmdResults)engine.Run();
            Assert.AreEqual(0, res.Mods.Count);

            ParentSpectrumMatch newPsm = new TestParentSpectrumMatch(588.22520189093 + 42);
            var parentProtein = new Protein("NNNNN", "accession", new Dictionary<int, List<MorpheusModification>>(), null, null, null, null, null, 0, false, false);
            var modPep = new PeptideWithPossibleModifications(1, 5, parentProtein, 0, "ugh");
            //var twoBasedVariableAndLocalizeableModificationss = new Dictionary<int, MorpheusModification>();
            List<MorpheusModification> variableModifications = new List<MorpheusModification>();
            var peptidesWithSetModifications = new HashSet<PeptideWithSetModifications> { modPep.GetPeptideWithSetModifications(variableModifications, 4096, 3).First() };
            Tolerance fragmentTolerance = null;
            IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile = null;
            var thisPSM = new PSMwithProteinHashSet(newPsm, peptidesWithSetModifications, fragmentTolerance, myMsDataFile, new List<ProductType> { ProductType.B, ProductType.Y });
            var newPsmWithFDR = new NewPsmWithFdr(thisPSM, 1, 0, 0);
            allResultingIdentifications.Add(newPsmWithFDR);

            engine = new GptmdEngine(allResultingIdentifications, isotopeErrors, gptmdModifications, combos, tol);
            res = (GptmdResults)engine.Run();
            Assert.AreEqual(1, res.Mods.Count);
            Assert.AreEqual(5, res.Mods["accession"].Count);
        }

        #endregion Public Methods

        #region Private Classes

        private class TestParentSpectrumMatch : ParentSpectrumMatch
        {

            #region Public Constructors

            public TestParentSpectrumMatch(double scanPrecursorMass) : base(null, double.NaN, double.NaN, scanPrecursorMass, 0, 0, 0, double.NaN, double.NaN, double.NaN)
            {
            }

            #endregion Public Constructors

            #region Public Methods

            public override CompactPeptide GetCompactPeptide(List<MorpheusModification> variableModifications, List<MorpheusModification> localizeableModifications)
            {
                throw new NotImplementedException();
            }

            #endregion Public Methods

        }

        #endregion Private Classes

    }
}