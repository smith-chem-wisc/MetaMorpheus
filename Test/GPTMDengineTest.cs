using InternalLogicEngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using OldInternalLogic;
using Spectra;
using System;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    public class GPTMDengineTest
    {
        #region Public Methods

        [Test]
        public void TestGPTMDengine()
        {
            List<NewPsmWithFDR> allResultingIdentifications = null;
            var gptmdModifications = new List<MorpheusModification> { new MorpheusModification("name", ModificationType.AminoAcidResidue, 'N', 42, null, null, '\0', double.NaN, false, null) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>>();
            double tol = 0.1;
            bool isotopeErrors = false;

            allResultingIdentifications = new List<NewPsmWithFDR>();
            var engine = new GPTMDEngine(allResultingIdentifications, isotopeErrors, gptmdModifications, combos, tol);
            var res = (GPTMDResults)engine.Run();
            Assert.AreEqual(0, res.mods.Count);

            ParentSpectrumMatch newPsm = new TestParentSpectrumMatch(588.22520189093 + 42);
            var parentProtein = new Protein("NNNNN", "accession", null, new Dictionary<int, List<MorpheusModification>>(), null, null, null, null, null, 0, false);
            var modPep = new PeptideWithPossibleModifications(1, 5, parentProtein, 0, "ugh");
            var twoBasedVariableAndLocalizeableModificationss = new Dictionary<int, MorpheusModification>();
            var peptidesWithSetModifications = new HashSet<PeptideWithSetModifications> { new PeptideWithSetModifications(modPep, twoBasedVariableAndLocalizeableModificationss) };
            Tolerance fragmentTolerance = null;
            IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile = null;
            var thisPSM = new PSMwithTargetDecoyKnown(newPsm, peptidesWithSetModifications, fragmentTolerance, myMsDataFile);
            var newPsmWithFDR = new NewPsmWithFDR(thisPSM, 1, 0, 0);
            allResultingIdentifications.Add(newPsmWithFDR);

            engine = new GPTMDEngine(allResultingIdentifications, isotopeErrors, gptmdModifications, combos, tol);
            res = (GPTMDResults)engine.Run();
            Assert.AreEqual(1, res.mods.Count);
            Assert.AreEqual(5, res.mods["accession"].Count);
        }

        #endregion Public Methods

        #region Private Classes

        private class TestParentSpectrumMatch : ParentSpectrumMatch
        {
            #region Public Constructors

            public TestParentSpectrumMatch(double scanPrecursorMass)
            {
                this.scanPrecursorMass = scanPrecursorMass;
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