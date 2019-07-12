using EngineLayer;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;

namespace Test
{
    [TestFixture]
    public static class SeqCoverageTest
    {
        [Test]
        public static void TryFailSequenceCoverage()
        {
            var prot1 = new Protein("MMKMMK", "prot1");

            ModificationMotif.TryGetMotif("M", out ModificationMotif motifM);
            Modification mod1 = new Modification(_originalId: "mod1", _modificationType: "mt", _target: motifM, _locationRestriction: "N-terminal.", _monoisotopicMass: 10);
            Modification mod2 = new Modification(_originalId: "mod2", _modificationType: "mt", _target: motifM, _locationRestriction: "Peptide N-terminal.", _monoisotopicMass: 10);
            Modification mod3 = new Modification(_originalId: "mod3", _modificationType: "mt", _target: motifM, _locationRestriction: "Anywhere.", _monoisotopicMass: 10);
            ModificationMotif.TryGetMotif("K", out ModificationMotif motifK);
            Modification mod4 = new Modification(_originalId: "mod4", _modificationType: "mt", _target: motifK, _locationRestriction: "Peptide C-terminal.", _monoisotopicMass: 10);
            Modification mod5 = new Modification(_originalId: "mod5", _modificationType: "mt", _target: motifK, _locationRestriction: "C-terminal.", _monoisotopicMass: 10);

            Dictionary<int, Modification> modsFor1 = new Dictionary<int, Modification>
            {
                {1, mod1},
                {3, mod3},
                {5, mod4},
            };
            Dictionary<int, Modification> modsFor2 = new Dictionary<int, Modification>
            {
                {1, mod2},
                {5, mod5},
            };
            Dictionary<int, Modification> modsFor3 = new Dictionary<int, Modification>
            {
                {1, mod1},
                {5, mod3},
                {8, mod5}
            };

            DigestionParams digestionParams = new DigestionParams();
            var pwsm1 = new PeptideWithSetModifications(prot1, digestionParams, 1, 3, CleavageSpecificity.Unknown, "",  0,  modsFor1,  0);
            var pwsm2 = new PeptideWithSetModifications(prot1, digestionParams, 4, 6, CleavageSpecificity.Unknown, "",  0,  modsFor2,  0);
            var pwsm3 = new PeptideWithSetModifications(prot1, digestionParams, 1, 6, CleavageSpecificity.Unknown, "",  0,  modsFor3,  0);

            HashSet<PeptideWithSetModifications> peptides = new HashSet<PeptideWithSetModifications>
            {
                pwsm1,
                pwsm2,
                pwsm3,
            };

            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""), 0, 0, "", new CommonParameters());

            var psm1 = new PeptideSpectralMatch(pwsm1, 0, 1, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            psm1.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);

            var psm2 = new PeptideSpectralMatch(pwsm2, 0, 1, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            psm2.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);

            var psm3 = new PeptideSpectralMatch(pwsm3, 0, 1, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            psm3.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);


            List<PeptideSpectralMatch> newPsms = new List<PeptideSpectralMatch>
            {
                psm1,
                psm2,
                psm3,
            };

            newPsms.ForEach(p => p.ResolveAllAmbiguities());

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(newPsms, true, new CommonParameters(), new List<string>());
            ProteinParsimonyResults fjkd = (ProteinParsimonyResults)ppe.Run();

            ProteinScoringAndFdrEngine psafe = new ProteinScoringAndFdrEngine(fjkd.ProteinGroups, newPsms, true, true, true, new CommonParameters(), new List<string>());

            psafe.Run();

            fjkd.ProteinGroups.First().CalculateSequenceCoverage();

            var firstSequenceCoverageDisplayList = fjkd.ProteinGroups.First().SequenceCoverageDisplayList.First();
            Assert.AreEqual("MMKMMK", firstSequenceCoverageDisplayList);
            var firstSequenceCoverageDisplayListWithMods = fjkd.ProteinGroups.First().SequenceCoverageDisplayListWithMods.First();
            Assert.AreEqual("[mod1 on M]-MM[mod3 on M]KM[mod3 on M]MK-[mod5 on K]", firstSequenceCoverageDisplayListWithMods);

            var firstModInfo = fjkd.ProteinGroups.First().ModsInfo.First();
            Assert.IsTrue(firstModInfo.Contains(@"#aa1[mod1 on M,info:occupancy=1.00(2/2)]"));
            Assert.IsTrue(firstModInfo.Contains(@"#aa2[mod3 on M,info:occupancy=0.50(1/2)]"));
            Assert.IsFalse(firstModInfo.Contains(@"#aa3"));
            Assert.IsTrue(firstModInfo.Contains(@"#aa4[mod3 on M,info:occupancy=0.50(1/2)]"));
            Assert.IsFalse(firstModInfo.Contains(@"#aa5"));
            Assert.IsTrue(firstModInfo.Contains(@"#aa6[mod5 on K,info:occupancy=1.00(2/2)]"));
        }
    }
}