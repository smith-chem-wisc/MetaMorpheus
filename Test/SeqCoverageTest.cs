using EngineLayer;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test
{
    [TestFixture]
    public static class SeqCoverageTest
    {
        #region Public Methods

        [Test]
        public static void TryFailSequenceCoverage()
        {
            var prot1 = new Protein("MMKMMK", "prot1");

            ModificationMotif.TryGetMotif("M", out ModificationMotif motifM);
            ModificationWithMass mod1 = new ModificationWithMass("mod1", "mt", motifM, TerminusLocalization.NProt, 10);
            ModificationWithMass mod2 = new ModificationWithMass("mod2", "mt", motifM, TerminusLocalization.NPep, 10);
            ModificationWithMass mod3 = new ModificationWithMass("mod3", "mt", motifM, TerminusLocalization.Any, 10);
            ModificationMotif.TryGetMotif("K", out ModificationMotif motifK);
            ModificationWithMass mod4 = new ModificationWithMass("mod4", "mt", motifK, TerminusLocalization.PepC, 10);
            ModificationWithMass mod5 = new ModificationWithMass("mod5", "mt", motifK, TerminusLocalization.ProtC, 10);

            Dictionary<int, ModificationWithMass> modsFor1 = new Dictionary<int, ModificationWithMass>
            {
                {1, mod1},
                {3, mod3},
                {5, mod4},
            };
            Dictionary<int, ModificationWithMass> modsFor2 = new Dictionary<int, ModificationWithMass>
            {
                {1, mod2},
                {5, mod5},
            };
            Dictionary<int, ModificationWithMass> modsFor3 = new Dictionary<int, ModificationWithMass>
            {
                {1, mod1},
                {5, mod3},
                {8, mod5}
            };

            var pwsm1 = new PeptideWithSetModifications(0, prot1, 1, 3, modsFor1);
            var pwsm2 = new PeptideWithSetModifications(0, prot1, 4, 6, modsFor2);
            var pwsm3 = new PeptideWithSetModifications(0, prot1, 1, 6, modsFor3);
            HashSet<PeptideWithSetModifications> peptides = new HashSet<PeptideWithSetModifications>
            {
                pwsm1,
                pwsm2,
                pwsm3,
            };

            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>
            {
                { pwsm1.CompactPeptide(TerminusType.None), new HashSet<PeptideWithSetModifications>{ pwsm1 } },
                { pwsm2.CompactPeptide(TerminusType.None), new HashSet<PeptideWithSetModifications>{ pwsm2 } },
                { pwsm3.CompactPeptide(TerminusType.None), new HashSet<PeptideWithSetModifications>{ pwsm3 } },
            };
            var digestionParams = new DigestionParams();
            IScan scan = new ThisTestScan();
            var psm1 = new PeptideSpectralMatch(pwsm1.CompactPeptide(TerminusType.None), 0, 1, 0, scan, digestionParams);
            psm1.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
            psm1.MatchToProteinLinkedPeptides(matching);
            var psm2 = new PeptideSpectralMatch(pwsm2.CompactPeptide(TerminusType.None), 0, 1, 0, scan, digestionParams);
            psm2.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
            psm2.MatchToProteinLinkedPeptides(matching);
            var psm3 = new PeptideSpectralMatch(pwsm3.CompactPeptide(TerminusType.None), 0, 1, 0, scan, digestionParams);
            psm3.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
            psm3.MatchToProteinLinkedPeptides(matching);

            List<PeptideSpectralMatch> newPsms = new List<PeptideSpectralMatch>
            {
                psm1,
                psm2,
                psm3,
            };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(matching, new HashSet<DigestionParams> { digestionParams }, true, new List<string>());
            ProteinParsimonyResults fjkd = (ProteinParsimonyResults)ppe.Run();

            ProteinScoringAndFdrEngine psafe = new ProteinScoringAndFdrEngine(fjkd.ProteinGroups, newPsms, true, true, true, new List<string>());

            psafe.Run();

            fjkd.ProteinGroups.First().CalculateSequenceCoverage();

            var firstSequenceCoverageDisplayList = fjkd.ProteinGroups.First().SequenceCoverageDisplayList.First();
            Assert.AreEqual("MMKMMK", firstSequenceCoverageDisplayList);
            var firstSequenceCoverageDisplayListWithMods = fjkd.ProteinGroups.First().SequenceCoverageDisplayListWithMods.First();
            Assert.AreEqual("[mod1]-MM[mod3]KM[mod3]MK-[mod5]", firstSequenceCoverageDisplayListWithMods);

            var firstModInfo = fjkd.ProteinGroups.First().ModsInfo.First();
            Assert.IsTrue(firstModInfo.Contains(@"#aa1[mod1,info:occupancy=1.00(2/2)]"));
            Assert.IsTrue(firstModInfo.Contains(@"#aa2[mod3,info:occupancy=0.50(1/2)]"));
            Assert.IsFalse(firstModInfo.Contains(@"#aa3"));
            Assert.IsTrue(firstModInfo.Contains(@"#aa4[mod3,info:occupancy=0.50(1/2)]"));
            Assert.IsFalse(firstModInfo.Contains(@"#aa5"));
            Assert.IsTrue(firstModInfo.Contains(@"#aa6[mod5,info:occupancy=1.00(2/2)]"));
        }

        [Test]
        public static void MultipleProteaseSelectionTest()
        {
            Protein ParentProtein = new Protein("MOAT", "accession1");

            var protease = new Protease("TestProtease1", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("O", TerminusType.C), new Tuple<string, TerminusType>("T", TerminusType.N) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            GlobalVariables.ProteaseDictionary.Add(protease.Name, protease);
            DigestionParams multiProtease = new DigestionParams(protease: protease.Name, MaxMissedCleavages: 0, MinPeptideLength: 1, InitiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var digestedList = ParentProtein.Digest(multiProtease, new List<ModificationWithMass>(), new List<ModificationWithMass>()).ToList();
            var sequences = digestedList.Select(p => p.BaseSequence).ToList();
            Assert.That(sequences.Count == 3);
            Assert.That(sequences.Contains("MO"));
            Assert.That(sequences.Contains("A"));
            Assert.That(sequences.Contains("T"));
        }

        [Test]
        public static void MultipleProteaseSelectionTestMissedCleavage()
        {
            Protein ParentProtein = new Protein("MOAT", "accession1");

            var protease = new Protease("TestProtease2", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("O", TerminusType.C), new Tuple<string, TerminusType>("T", TerminusType.N) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            GlobalVariables.ProteaseDictionary.Add(protease.Name, protease);
            DigestionParams multiProtease = new DigestionParams(protease: protease.Name, MaxMissedCleavages: 1, MinPeptideLength: 1, InitiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var digestedList = ParentProtein.Digest(multiProtease, new List<ModificationWithMass>(), new List<ModificationWithMass>()).ToList();
            var sequences = digestedList.Select(p => p.BaseSequence).ToList();
            Assert.That(sequences.Count == 5);
            Assert.That(sequences.Contains("MOA"));
            Assert.That(sequences.Contains("AT"));
            Assert.That(sequences.Contains("MO"));
            Assert.That(sequences.Contains("A"));
            Assert.That(sequences.Contains("T"));
        }

        [Test]
        public static void MultipleProteaseSelectionTestPreventCleavage()
        {
            Protein ParentProtein = new Protein("MOAT", "accession1");

            var protease = new Protease("TestProtease3", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("O", TerminusType.C), new Tuple<string, TerminusType>("T", TerminusType.N) }, new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("A", TerminusType.C) }, CleavageSpecificity.Full, null, null, null);
            GlobalVariables.ProteaseDictionary.Add(protease.Name, protease);
            DigestionParams multiProtease = new DigestionParams(protease: protease.Name, MaxMissedCleavages: 0, MinPeptideLength: 1, InitiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var digestedList = ParentProtein.Digest(multiProtease, new List<ModificationWithMass>(), new List<ModificationWithMass>()).ToList();
            var sequences = digestedList.Select(p => p.BaseSequence).ToList();
            Assert.That(sequences.Count == 2);
            Assert.That(sequences.Contains("MOA"));
            Assert.That(sequences.Contains("T"));
        }

        [Test]
        public static void ReadCustomFile()
        {
            Protein ParentProtein = new Protein("OKAREDY", "accession1");
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DoubleProtease.tsv");
            Assert.That(File.Exists(path));

            var proteaseDict = GlobalVariables.LoadProteaseDictionary(path);
            Assert.That(proteaseDict.ContainsKey("Test1"));
            Assert.That(proteaseDict.ContainsKey("Test2"));
            Assert.That(proteaseDict.ContainsKey("Test3"));
            GlobalVariables.ProteaseDictionary.Add("Test1", proteaseDict["Test1"]);

            DigestionParams multiProtease1 = new DigestionParams(protease: "Test1", MaxMissedCleavages: 0, MinPeptideLength: 1, InitiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var digestedList1 = ParentProtein.Digest(multiProtease1, new List<ModificationWithMass>(), new List<ModificationWithMass>()).ToList();
            GlobalVariables.ProteaseDictionary.Remove("Test1");

            var sequences = digestedList1.Select(p => p.BaseSequence).ToList();
            Assert.That(sequences.Count == 3);
            Assert.That(sequences.Contains("OK"));
            Assert.That(sequences.Contains("A"));
            Assert.That(sequences.Contains("REDY"));

            GlobalVariables.ProteaseDictionary.Add("Test2", proteaseDict["Test2"]);
            DigestionParams multiProtease2 = new DigestionParams(protease: "Test2", MaxMissedCleavages: 0, MinPeptideLength: 1, InitiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var digestedList2 = ParentProtein.Digest(multiProtease2, new List<ModificationWithMass>(), new List<ModificationWithMass>()).ToList();
            GlobalVariables.ProteaseDictionary.Remove("Test2");
            var sequences2 = digestedList2.Select(p => p.BaseSequence).ToList();
            Assert.That(sequences2.Count == 3);
            Assert.That(sequences2.Contains("OK"));
            Assert.That(sequences2.Contains("ARED"));
            Assert.That(sequences2.Contains("Y"));

            GlobalVariables.ProteaseDictionary.Add("Test3", proteaseDict["Test3"]);
            DigestionParams multiProtease3 = new DigestionParams(protease: "Test3", MaxMissedCleavages: 0, MinPeptideLength: 1, InitiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var digestedList3 = ParentProtein.Digest(multiProtease3, new List<ModificationWithMass>(), new List<ModificationWithMass>()).ToList();
            GlobalVariables.ProteaseDictionary.Remove("Test3");
            var sequences3 = digestedList3.Select(p => p.BaseSequence).ToList();
            Assert.That(sequences3.Count == 2);
            Assert.That(sequences3.Contains("OK"));
            Assert.That(sequences3.Contains("AREDY"));
        }
        #endregion Public Methods
    }
}