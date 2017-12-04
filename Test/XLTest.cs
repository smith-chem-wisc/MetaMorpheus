using EngineLayer;
using EngineLayer.CrosslinkSearch;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;
using System;

namespace Test
{
    [TestFixture]
    public static class XLTest
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
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();

            var ye = prot.Digest(digestionParams, new List<ModificationWithMass>(), variableModifications).ToList();

            var pep = ye[0];
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

            var pep2 = ye[2];
            Assert.AreEqual("MNNNKQQQQ", pep2.BaseSequence);
            var n2 = pep2.CompactPeptide(TerminusType.None).NTerminalMasses;
            var c2 = pep2.CompactPeptide(TerminusType.None).CTerminalMasses;
            Assert.AreEqual(n2.Count(), 8);
            Assert.AreEqual(c2.Count(), 8);
            Assert.AreEqual(n2[4] - n2[3], 128.09496301518999, 1e-6);
            var x2 = PsmCross.xlPosCal(pep2.CompactPeptide(TerminusType.None), crosslinker).ToArray();
            Assert.AreEqual(x2[0], 4);
        }

        [Test]
        public static void XLGetRankArrayOfDoubleArray()
        {
            double[] mz = new double[] { 1.0, 1.3, 1.5, 1.7, 1.9, 2.1 };
            double[] intensity = new double[] { 1.1, 1.1, 0.5, 3.2, 0.5, 6.0};
            int[] rank = PsmCross.GenerateIntensityRanks(mz, intensity);
            int[] Rank = new int[] { 4, 3, 6, 2, 5, 1 };
            Assert.AreEqual(rank, Rank);
        }

        #endregion Public Methods
    }
}