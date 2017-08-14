﻿using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class TestPsm
    {
        #region Public Methods

        [Test]
        public static void TestPsmHeader()
        {
            PeptideWithSetModifications pepWithSetMods = new Protein("MQQQQQQQ", "accession1").Digest(new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null), 0, null, null, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>()).First().GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).First();
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(pepWithSetMods, "quadratic");
            IMzPeak peak = new MzPeak(4, 4);
            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> scann = myMsDataFile.GetOneBasedScan(2) as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>;
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(scann, peak, 1, null);
            Psm psm = new Psm(pepWithSetMods.CompactPeptide(TerminusType.None), 1, 2, 3, scan);

            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), Psm.GetTabSeparatedHeader().Count(f => f == '\t'));

            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>
            {
                { pepWithSetMods.CompactPeptide(TerminusType.None), new HashSet<PeptideWithSetModifications>{ pepWithSetMods } }
            };

            psm.MatchToProteinLinkedPeptides(matching);

            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), Psm.GetTabSeparatedHeader().Count(f => f == '\t'));

            Tolerance fragmentTolerance = new PpmTolerance(10);
            List<ProductType> lp = new List<ProductType> { ProductType.B };

            new LocalizationEngine(new List<Psm> { psm }, lp, myMsDataFile, fragmentTolerance, new List<string>(), false).Run();

            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), Psm.GetTabSeparatedHeader().Count(f => f == '\t'));

            psm.SetFdrValues(6, 6, 6, 6, 6, 6);

            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), Psm.GetTabSeparatedHeader().Count(f => f == '\t'));
        }

        #endregion Public Methods
    }
}