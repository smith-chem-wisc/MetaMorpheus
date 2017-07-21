using Chemistry;
using EngineLayer;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;

using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;
using static Chemistry.PeriodicTable;
using EngineLayer.CrosslinkSearch;


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
            var ye = prot.Digest(protease, 1, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).ToList();

            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();

            var pep = ye[0].GetPeptidesWithSetModifications(variableModifications, 4096, 3).First();
            Assert.AreEqual(pep.BaseSequence, "MNNNK");
            CrosslinkerTypeClass crosslinker = new CrosslinkerTypeClass();
            crosslinker.SelectCrosslinker(CrosslinkerType.DSS);
            Assert.AreEqual(crosslinker.CrosslinkerModSite, 'K');
            Assert.AreEqual(Residue.GetResidue(crosslinker.CrosslinkerModSite).MonoisotopicMass, 128.09496301518999, 1e-9);
            var n = pep.CompactPeptide.NTerminalMasses;
            var c = pep.CompactPeptide.CTerminalMasses;
            Assert.AreEqual(n.Count(), 4);
            Assert.AreEqual(c.Count(), 4);
            Assert.AreEqual(c[0], 128.09496301518999);
            var x = PsmCross.xlPosCal(pep.CompactPeptide, crosslinker).ToArray();
            Assert.AreEqual(x[0], 4);

            var pep2 = ye[2].GetPeptidesWithSetModifications(variableModifications, 4096, 3).First();
            Assert.AreEqual(pep2.BaseSequence, "MNNNKQQQQ");
            var n2 = pep2.CompactPeptide.NTerminalMasses;
            var c2 = pep2.CompactPeptide.CTerminalMasses;
            Assert.AreEqual(n2.Count(), 8);
            Assert.AreEqual(c2.Count(), 8);
            Assert.AreEqual(n2[4] - n2[3] , 128.09496301518999, 1e-9);
            var x2 = PsmCross.xlPosCal(pep2.CompactPeptide, crosslinker).ToArray();
            Assert.AreEqual(x2[0], 4);
        }

        #endregion
    }
}
