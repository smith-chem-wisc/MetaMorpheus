﻿using EngineLayer;
using EngineLayer.Localization;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class TestPsm
    {
        #region Public Methods

        [Test]
        public static void TestPsmHeader()
        {
            DigestionParams digestionParams = new DigestionParams();
            PeptideWithSetModifications pepWithSetMods = new Protein("MQQQQQQQ", "accession1").Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).First();
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(pepWithSetMods, "quadratic");
            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> scann = myMsDataFile.GetOneBasedScan(2) as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>;
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(scann, 4, 1, null);
            PeptideSpectralMatch psm = new PeptideSpectralMatch(pepWithSetMods.CompactPeptide(TerminusType.None), 1, 2, 3, scan);

            var t = psm.ToString();
            var tabsepheader = PeptideSpectralMatch.GetTabSeparatedHeader();
            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));

            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>
            {
                { pepWithSetMods.CompactPeptide(TerminusType.None), new HashSet<PeptideWithSetModifications>{ pepWithSetMods } }
            };

            psm.MatchToProteinLinkedPeptides(matching);

            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));

            Tolerance fragmentTolerance = new PpmTolerance(10);
            List<ProductType> lp = new List<ProductType> { ProductType.B };

            new LocalizationEngine(new List<PeptideSpectralMatch> { psm }, lp, myMsDataFile, fragmentTolerance, new List<string>(), false).Run();

            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));

            psm.SetFdrValues(6, 6, 6, 6, 6, 6, 0, 0, 0, false);

            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));
        }

        #endregion Public Methods
    }
}