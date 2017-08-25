using EngineLayer;
using EngineLayer.CrosslinkSearch;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class XLTest
    {
        #region Public Methods

        [Test]
        public static void XLTestXlPosCal()
        {
            var prot = new Protein("MNNNKQQQQ", null);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            DigestionParams digestionParams = new DigestionParams
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain,
                MaxMissedCleavages = 2,
                Protease = protease,
                MinPeptideLength = 1
            };
            var ye = prot.Digest(digestionParams, new List<ModificationWithMass>()).ToList();

            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();

            var pep = ye[0].GetPeptidesWithSetModifications(digestionParams, variableModifications).First();
            Assert.AreEqual(pep.BaseSequence, "MNNNK");
            CrosslinkerTypeClass crosslinker = new CrosslinkerTypeClass();
            crosslinker.SelectCrosslinker(CrosslinkerType.DSS);
            Assert.AreEqual(crosslinker.CrosslinkerModSite, 'K');
            Assert.AreEqual(Residue.GetResidue(crosslinker.CrosslinkerModSite).MonoisotopicMass, 128.09496301518999, 1e-9);
            var n = pep.CompactPeptide(TerminusType.None).NTerminalMasses;
            var c = pep.CompactPeptide(TerminusType.None).CTerminalMasses;
            Assert.AreEqual(n.Count(), 4);
            Assert.AreEqual(c.Count(), 4);
            Assert.AreEqual(c[0], 128.09496301518999, 1e-6);
            var x = PsmCross.xlPosCal(pep.CompactPeptide(TerminusType.None), crosslinker).ToArray();
            Assert.AreEqual(x[0], 4);

            var pep2 = ye[2].GetPeptidesWithSetModifications(digestionParams, variableModifications).First();
            Assert.AreEqual("MNNNKQQQQ", pep2.BaseSequence);
            var n2 = pep2.CompactPeptide(TerminusType.None).NTerminalMasses;
            var c2 = pep2.CompactPeptide(TerminusType.None).CTerminalMasses;
            Assert.AreEqual(n2.Count(), 8);
            Assert.AreEqual(c2.Count(), 8);
            Assert.AreEqual(n2[4] - n2[3], 128.09496301518999, 1e-6);
            var x2 = PsmCross.xlPosCal(pep2.CompactPeptide(TerminusType.None), crosslinker).ToArray();
            Assert.AreEqual(x2[0], 4);
        }

        #endregion Public Methods
    }
}