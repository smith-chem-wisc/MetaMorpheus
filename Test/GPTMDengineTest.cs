using EngineLayer;
using EngineLayer.Gptmd;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
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
            ModificationMotif motifN;
            ModificationMotif.TryGetMotif("N", out motifN);
            var gptmdModifications = new List<ModificationWithMass> { new ModificationWithMass("21", null, motifN, ModificationSites.Any, 21.981943, null, new List<double> { 0 }, new List<double> { 21.981943 }, null) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>>();
            Tolerance precursorMassTolerance = new Tolerance(ToleranceUnit.PPM, 10);
            bool isotopeErrors = false;

            allResultingIdentifications = new List<NewPsmWithFdr>();
            var engine = new GptmdEngine(allResultingIdentifications, isotopeErrors, gptmdModifications, combos, precursorMassTolerance);
            var res = (GptmdResults)engine.Run();
            Assert.AreEqual(0, res.Mods.Count);

            PsmParent newPsm = new TestParentSpectrumMatch(588.22520189093 + 21.981943);
            var parentProtein = new Protein("NNNNN", "accession", null, new Dictionary<int, List<Modification>>(), null, null, null, null, null, false, false, null);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var modPep = parentProtein.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>()).First();
            //var twoBasedVariableAndLocalizeableModificationss = new Dictionary<int, MorpheusModification>();
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            var peptidesWithSetModifications = new HashSet<PeptideWithSetModifications> { modPep.GetPeptidesWithSetModifications(variableModifications, 4096, 3).First() };
            Tolerance fragmentTolerance = null;
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = null;
            var thisPSM = new PsmWithMultiplePossiblePeptides(newPsm, peptidesWithSetModifications, fragmentTolerance, myMsDataFile, new List<ProductType> { ProductType.B, ProductType.Y });

            NewPsmWithFdr thePsmwithfdr = new NewPsmWithFdr(thisPSM);
            thePsmwithfdr.SetValues(1, 0, 0, 1, 0, 0);
            allResultingIdentifications.Add(thePsmwithfdr);

            engine = new GptmdEngine(allResultingIdentifications, isotopeErrors, gptmdModifications, combos, precursorMassTolerance);
            res = (GptmdResults)engine.Run();
            Assert.AreEqual(1, res.Mods.Count);
            Assert.AreEqual(5, res.Mods["accession"].Count);
        }

        [Test]
        public static void TestCombos()
        {
            List<NewPsmWithFdr> allIdentifications = null;
            ModificationMotif motifN;
            ModificationMotif.TryGetMotif("N", out motifN);
            ModificationMotif motifP;
            ModificationMotif.TryGetMotif("P", out motifP);
            var gptmdModifications = new List<ModificationWithMass> { new ModificationWithMass("21", null, motifN, ModificationSites.Any, 21.981943,null, new List<double> { 0 }, new List<double> { 21.981943 },  null),
                                                                      new ModificationWithMass("16", null, motifP, ModificationSites.Any, 15.994915,null, new List<double> { 0 }, new List<double> { 15.994915 },  null) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>> { new Tuple<double, double>(21.981943, 15.994915) };
            Tolerance precursorMassTolerance = new Tolerance(ToleranceUnit.PPM, 10);
            bool isotopeErrors = false;
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            PsmParent newPsm = new TestParentSpectrumMatch(651.297638557 + 21.981943 + 15.994915);
            var parentProtein = new Protein("NNNPPP", "accession", null, new Dictionary<int, List<Modification>>(), null, null, null, null, null, false, false, null);
            var modPep = parentProtein.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>()).First();

            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            var peptidesWithSetModifications = new HashSet<PeptideWithSetModifications> { modPep.GetPeptidesWithSetModifications(variableModifications, 4096, 3).First() };
            Tolerance fragmentTolerance = null;
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = null;
            var thisPSM = new PsmWithMultiplePossiblePeptides(newPsm, peptidesWithSetModifications, fragmentTolerance, myMsDataFile, new List<ProductType> { ProductType.B, ProductType.Y });
            NewPsmWithFdr thePsmwithfdr = new NewPsmWithFdr(thisPSM);
            thePsmwithfdr.SetValues(1, 0, 0, 1, 0, 0);
            allIdentifications = new List<NewPsmWithFdr> { thePsmwithfdr };

            var engine = new GptmdEngine(allIdentifications, isotopeErrors, gptmdModifications, combos, precursorMassTolerance);
            var res = (GptmdResults)engine.Run();
            Assert.AreEqual(1, res.Mods.Count);
            Assert.AreEqual(6, res.Mods["accession"].Count);
            Assert.AreEqual(3, res.Mods["accession"].Where(b => b.Item2.id.Equals("21")).Count());
            Assert.AreEqual(3, res.Mods["accession"].Where(b => b.Item2.id.Equals("16")).Count());
        }

        #endregion Public Methods

        #region Private Classes

        private class TestParentSpectrumMatch : PsmParent
        {

            #region Public Constructors

            public TestParentSpectrumMatch(double scanPrecursorMass) : base(null, double.NaN, double.NaN, scanPrecursorMass, 0, 0, 0, 0, double.NaN, double.NaN, double.NaN, 1)
            {
            }

            #endregion Public Constructors

            #region Public Methods

            public override CompactPeptide GetCompactPeptide(List<ModificationWithMass> variableModifications, List<ModificationWithMass> localizeableModifications, List<ModificationWithMass> fixedModifications)
            {
                throw new NotImplementedException();
            }

            #endregion Public Methods

        }

        #endregion Private Classes

    }
}