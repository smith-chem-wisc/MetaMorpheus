using System;
using EngineLayer;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using FlashLFQ;
using static iText.IO.Image.Jpeg2000ImageData;
using TaskLayer;

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

            var psm1 = new PeptideSpectralMatch(pwsm1, 0, 1, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            psm1.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);

            var psm2 = new PeptideSpectralMatch(pwsm2, 0, 1, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            psm2.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);

            var psm3 = new PeptideSpectralMatch(pwsm3, 0, 1, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            psm3.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);


            List<PeptideSpectralMatch> newPsms = new List<PeptideSpectralMatch>
            {
                psm1,
                psm2,
                psm3,
            };

            newPsms.ForEach(p => p.ResolveAllAmbiguities());

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(newPsms, true, new CommonParameters(), null, new List<string>());
            ProteinParsimonyResults fjkd = (ProteinParsimonyResults)ppe.Run();

            ProteinScoringAndFdrEngine psafe = new ProteinScoringAndFdrEngine(fjkd.ProteinGroups, newPsms, true, true, true, new CommonParameters(), null, new List<string>());

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
            Console.WriteLine("Test output: " + firstSequenceCoverageDisplayList);
        }


        [Test]
        public static void TestFragmentSequenceCoverage()
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
            var pwsm1 = new PeptideWithSetModifications(prot1, digestionParams, 1, 3, CleavageSpecificity.Unknown, "", 0, modsFor1, 0);
            var pwsm2 = new PeptideWithSetModifications(prot1, digestionParams, 4, 6, CleavageSpecificity.Unknown, "", 0, modsFor2, 0);
            var pwsm3 = new PeptideWithSetModifications(prot1, digestionParams, 1, 6, CleavageSpecificity.Unknown, "", 0, modsFor3, 0);



            Product productb1 = new Product(ProductType.b, 0, 0, 1, 1, 0);
            Product productb2 = new Product(ProductType.b, 0, 0, 2, 2, 0);
            Product productb4 = new Product(ProductType.b, 0, 0, 4, 4, 0);
            Product producty1 = new Product(ProductType.y, 0, 0, 1, 0, 0);
            Product producty2 = new Product(ProductType.y, 0, 0, 2, 0, 0);

            MatchedFragmentIon mfib1 = new MatchedFragmentIon(ref productb1, 0, 0, 1);
            MatchedFragmentIon mfib2 = new MatchedFragmentIon(ref productb2, 0, 0, 2);
            MatchedFragmentIon mfib4 = new MatchedFragmentIon(ref productb4, 0, 0, 2);
            MatchedFragmentIon mfiy1 = new MatchedFragmentIon(ref producty1, 0, 0, 2);
            MatchedFragmentIon mfiy2 = new MatchedFragmentIon(ref producty2, 0, 0, 2);

            List<MatchedFragmentIon> mfis1 = new List<MatchedFragmentIon> { mfib1, mfib2};

            List<MatchedFragmentIon> mfis2 = new List<MatchedFragmentIon> { mfib1, mfib2, mfib4, mfiy2 };

            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""), 0, 0, "", new CommonParameters());

            var psm1 = new PeptideSpectralMatch(pwsm1, 0, 1, 0, scan, new CommonParameters(), mfis1);
            psm1.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);

            var psm2 = new PeptideSpectralMatch(pwsm2, 0, 1, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            psm2.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);

            var psm3 = new PeptideSpectralMatch(pwsm3, 0, 1, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            psm3.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);


            List<PeptideSpectralMatch> newPsms = new List<PeptideSpectralMatch>
            {
                psm1,
                psm2,
                psm3,
            };

            newPsms.ForEach(p => p.ResolveAllAmbiguities());

            foreach (var psm in newPsms)
            {
                     psm.AddFragmentCoveragePSMs();
            }

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(newPsms, true, new CommonParameters(), null, new List<string>());
            ProteinParsimonyResults fjkd = (ProteinParsimonyResults)ppe.Run();

            ProteinScoringAndFdrEngine psafe = new ProteinScoringAndFdrEngine(fjkd.ProteinGroups, newPsms, true, true, true, new CommonParameters(), null, new List<string>());

            psafe.Run();

            foreach (var group in fjkd.ProteinGroups)
            {
                group.CalculateSequenceCoverage();
            }

            // fjkd.ProteinGroups.First().CalculateSequenceCoverage();

            Console.WriteLine("-------");

            var firstSequenceCoverageDisplayList = fjkd.ProteinGroups.First().FragmentSequenceCoverageDisplayList.First();
            Console.WriteLine("Test output: " + firstSequenceCoverageDisplayList);

            


           // var prot1 = new Protein("PEPTIDEPEPTIDEKPEPTIDEPEPTIDEKPEPKPEPEPEPKTIDETIDERTIDEK", "prot1");

           // ModificationMotif.TryGetMotif("M", out ModificationMotif motifM);
           // Modification mod1 = new Modification(_originalId: "mod1", _modificationType: "mt", _target: motifM, _locationRestriction: "N-terminal.", _monoisotopicMass: 10);
           // Modification mod2 = new Modification(_originalId: "mod2", _modificationType: "mt", _target: motifM, _locationRestriction: "Peptide N-terminal.", _monoisotopicMass: 10);
           // Modification mod3 = new Modification(_originalId: "mod3", _modificationType: "mt", _target: motifM, _locationRestriction: "Anywhere.", _monoisotopicMass: 10);
           // ModificationMotif.TryGetMotif("K", out ModificationMotif motifK);
           // Modification mod4 = new Modification(_originalId: "mod4", _modificationType: "mt", _target: motifK, _locationRestriction: "Peptide C-terminal.", _monoisotopicMass: 10);
           // Modification mod5 = new Modification(_originalId: "mod5", _modificationType: "mt", _target: motifK, _locationRestriction: "C-terminal.", _monoisotopicMass: 10);

           // Dictionary<int, Modification> modsFor1 = new Dictionary<int, Modification>
           // {
           //     {1, mod1},
           //     {3, mod3},
           //     {5, mod4},
           // };
           // Dictionary<int, Modification> modsFor2 = new Dictionary<int, Modification>
           // {
           //     {1, mod2},
           //     {5, mod5},
           // };
           // Dictionary<int, Modification> modsFor3 = new Dictionary<int, Modification>
           // {
           //     {1, mod1},
           //     {5, mod3},
           //     {8, mod5}
           // };

           // DigestionParams digestionParams = new DigestionParams();
           // var pwsm1 = new PeptideWithSetModifications(prot1, digestionParams, 1, 15, CleavageSpecificity.Unknown, "", 0, modsFor1, 0);
           // var pwsm2 = new PeptideWithSetModifications(prot1, digestionParams, 1, 30, CleavageSpecificity.Unknown, "", 1, modsFor2, 0);
           // var pwsm3 = new PeptideWithSetModifications(prot1, digestionParams, 1, 15, CleavageSpecificity.Unknown, "", 0, modsFor3, 0);
           // var pwsm4 = new PeptideWithSetModifications(prot1, digestionParams, 31, 42, CleavageSpecificity.Unknown, "", 0, modsFor3, 0);

           // Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
           //     0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""), 0, 0, "", new CommonParameters());

           // //Protein p1 = new Protein("PEPTIDEPEPTIDEKPEPTIDEPEPTIDEKPEPKPEPEPEPKTIDETIDERTIDEK", null);
           // //CommonParameters commonParameters = new CommonParameters();
           // //PeptideWithSetModifications pep1 = p1.Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()).ToList()[0];

           // //Protein p2 = new Protein("PEPTIDEPEPTIDEKPEPTIDEPEPTIDEKPEPKPEPEPEPKTIDETIDERTIDEK", null);
           // //PeptideWithSetModifications pep2 = p2.Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()).ToList()[4];

           // //Protein p3 = new Protein("PEPTIDEPEPTIDEKPEPTIDEPEPTIDEKPEPKPEPEPEPKTIDETIDERTIDEK", null);
           // //PeptideWithSetModifications pep3 = p3.Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()).ToList()[0];

           // //Protein p4 = new Protein("PEPTIDEPEPTIDEKPEPTIDEPEPTIDEKPEPKPEPEPEPKTIDETIDERTIDEK", null);
           // //PeptideWithSetModifications pep4 = p4.Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()).ToList()[6];

           // // TestDataFile t = new TestDataFile(new List<PeptideWithSetModifications> { pep1, pep2, pep3, pep4 });

           // Product productb1 = new Product(ProductType.b, 0, 0, 1, 1, 0);
           // Product productb2 = new Product(ProductType.b, 0, 0, 2, 2, 0);
           // Product productb3 = new Product(ProductType.b, 0, 0, 3, 3, 0);
           // Product productb4 = new Product(ProductType.b, 0, 0, 4, 4, 0);
           // Product productb5 = new Product(ProductType.b, 0, 0, 5, 5, 0);
           // Product productb6 = new Product(ProductType.b, 0, 0, 6, 6, 0);


           // Product producty1 = new Product(ProductType.y, 0, 0, 1, 0, 0);
           // Product producty2 = new Product(ProductType.y, 0, 0, 2, 0, 0);
           // Product producty3 = new Product(ProductType.y, 0, 0, 3, 0, 0);
           // Product producty4 = new Product(ProductType.y, 0, 0, 4, 0, 0);
           // Product producty5 = new Product(ProductType.y, 0, 0, 5, 0, 0);
           // Product producty6 = new Product(ProductType.y, 0, 0, 6, 0, 0);
           // Product producty7 = new Product(ProductType.y, 0, 0, 7, 0, 0);


           // MatchedFragmentIon mfib1 = new MatchedFragmentIon(ref productb1, 0, 0, 1);
           // MatchedFragmentIon mfib2 = new MatchedFragmentIon(ref productb2, 0, 0, 2);
           // MatchedFragmentIon mfib3 = new MatchedFragmentIon(ref productb3, 0, 0, 3);
           // MatchedFragmentIon mfib4 = new MatchedFragmentIon(ref productb4, 0, 0, 1);
           // MatchedFragmentIon mfib5 = new MatchedFragmentIon(ref productb5, 0, 0, 1);
           // MatchedFragmentIon mfib6 = new MatchedFragmentIon(ref productb6, 0, 0, 1);

           // MatchedFragmentIon mfiy1 = new MatchedFragmentIon(ref producty1, 0, 0, 2);
           // MatchedFragmentIon mfiy2 = new MatchedFragmentIon(ref producty2, 0, 0, 2);
           // MatchedFragmentIon mfiy3 = new MatchedFragmentIon(ref producty3, 0, 0, 1);
           // MatchedFragmentIon mfiy4 = new MatchedFragmentIon(ref producty4, 0, 0, 1);
           // MatchedFragmentIon mfiy5 = new MatchedFragmentIon(ref producty5, 0, 0, 1);
           // MatchedFragmentIon mfiy6 = new MatchedFragmentIon(ref producty6, 0, 0, 1);
           // MatchedFragmentIon mfiy7 = new MatchedFragmentIon(ref producty7, 0, 0, 1);

           // List<MatchedFragmentIon> mfis1 = new List<MatchedFragmentIon> { mfib1, mfib2, mfib4, mfib5, mfib6, mfiy1, mfiy2, mfiy3, mfiy4, mfiy5, mfiy6 };
           // List<MatchedFragmentIon> mfis2 = new List<MatchedFragmentIon> { mfib1, mfib2, mfib4, mfib5, mfib6, mfiy1, mfiy2, mfiy3, mfiy4, mfiy5, mfiy6, mfiy7 };
           // List<MatchedFragmentIon> mfis3 = new List<MatchedFragmentIon> { mfib1, mfib2, mfib3, mfib4, mfib5, mfib6, mfiy1, mfiy2, mfiy3, mfiy4, mfiy5, mfiy6 };
           // List<MatchedFragmentIon> mfis4 = new List<MatchedFragmentIon> { mfib1, mfib2, mfib4, mfib5, mfib6, mfiy1, mfiy2, mfiy3, mfiy4, mfiy5, mfiy6, mfiy7 };

           // //MsDataScan mzLibScan1 = t.GetOneBasedScan(2);
           //// Ms2ScanWithSpecificMass scan1 = new Ms2ScanWithSpecificMass(mzLibScan1, 0, 1, null, new CommonParameters());
           // PeptideSpectralMatch psm1 = new PeptideSpectralMatch(pwsm1, 0, 0, 0, scan, new CommonParameters(), mfis1);

           // //MsDataScan mzLibScan2 = t.GetOneBasedScan(4);
           // //Ms2ScanWithSpecificMass scan2 = new Ms2ScanWithSpecificMass(mzLibScan2, 0, 1, null, new CommonParameters());
           // PeptideSpectralMatch psm2 = new PeptideSpectralMatch(pwsm2, 0, 0, 0, scan, new CommonParameters(), mfis2);

           // //MsDataScan mzLibScan3 = t.GetOneBasedScan(6);
           // //Ms2ScanWithSpecificMass scan3 = new Ms2ScanWithSpecificMass(mzLibScan3, 0, 1, null, new CommonParameters());
           // PeptideSpectralMatch psm3 = new PeptideSpectralMatch(pwsm3, 0, 0, 0, scan, new CommonParameters(), mfis3);

           // //MsDataScan mzLibScan4 = t.GetOneBasedScan(8);
           // //Ms2ScanWithSpecificMass scan4 = new Ms2ScanWithSpecificMass(mzLibScan4, 0, 1, null, new CommonParameters());
           // PeptideSpectralMatch psm4 = new PeptideSpectralMatch(pwsm4, 0, 0, 0, scan, new CommonParameters(), mfis4);

           // psm1.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0); // valid psm
           // psm1.ResolveAllAmbiguities();

           // psm2.SetFdrValues(0, 0, 0.02, 0, 0, 0, 0, 0); // psm above fdr cutoff
           // psm2.ResolveAllAmbiguities();

           // psm3.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0); // ambiguous peptide

           // psm4.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0); // valid psm for internal peptide
           // psm4.ResolveAllAmbiguities();

           // var allPsms = new List<PeptideSpectralMatch> { psm1, psm2, psm3, psm4 };

           // foreach (var psm in allPsms)
           // {
           //     psm.AddFragmentCoveragePSMs();
           // }

           // ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(allPsms, true, new CommonParameters(), null, new List<string>());
           // ProteinParsimonyResults fjkd = (ProteinParsimonyResults)ppe.Run();

           // ProteinScoringAndFdrEngine psafe = new ProteinScoringAndFdrEngine(fjkd.ProteinGroups, allPsms, false, true, true, new CommonParameters(), null, new List<string>());

           // psafe.Run();

           // fjkd.ProteinGroups.First().CalculateSequenceCoverage();
           // var firstSequenceCoverageDisplayList = fjkd.ProteinGroups.First().FragmentSequenceCoverageDisplayList.First();



           // Console.WriteLine(firstSequenceCoverageDisplayList);
            
        }
    }
}