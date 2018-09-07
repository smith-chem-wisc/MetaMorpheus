using EngineLayer;
using EngineLayer.Localization;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public static class TestPsm
    {
        [Test]
        public static void TestPsmHeader()
        {
            DigestionParams digestionParams = new DigestionParams();
            PeptideWithSetModifications pepWithSetMods = new Protein(
                "MQQQQQQQ",
                "accession1",
                "org",
                new List<Tuple<string, string>> { new Tuple<string, string>("geneNameType", "geneName") },
                new Dictionary<int, List<Modification>> { { 2, new List<Modification> { new Modification("mod", "mod") } } },
                name: "name",
                fullName: "fullName",
                sequenceVariations: new List<SequenceVariation> { new SequenceVariation(2, "P", "Q", "changed this sequence") })
                    .Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();
            MsDataFile myMsDataFile = new TestDataFile(pepWithSetMods, "quadratic");
            MsDataScan scann = myMsDataFile.GetOneBasedScan(2);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(scann, 4, 1, null);

            var theoreticalIons = pepWithSetMods.Fragment(DissociationType.HCD, FragmentationTerminus.Both).ToList();
            var matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan.TheScan.MassSpectrum, theoreticalIons, new CommonParameters());
            PeptideSpectralMatch psm = new PeptideSpectralMatch(pepWithSetMods, 1, 2, 3, scan, digestionParams, matchedIons);
            psm.ResolveAllAmbiguities();

            var t = psm.ToString();
            var tabsepheader = PeptideSpectralMatch.GetTabSeparatedHeader();
            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));

            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));

            Tolerance fragmentTolerance = new PpmTolerance(10);
            List<ProductType> lp = new List<ProductType> { ProductType.b };

            new LocalizationEngine(new List<PeptideSpectralMatch> { psm }, lp, myMsDataFile, new CommonParameters(productMassTolerance: fragmentTolerance), new List<string>()).Run();

            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));

            psm.SetFdrValues(6, 6, 6, 6, 6, 6, 0, 0, 0, true);

            Assert.AreEqual(psm.ToString().Count(f => f == '\t'), PeptideSpectralMatch.GetTabSeparatedHeader().Count(f => f == '\t'));
        }
    }
}