using Chemistry;
using EngineLayer;
using EngineLayer.CrosslinkSearch;
using EngineLayer.Indexing;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;
using MzLibUtil;
using Nett;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using OxyPlot.Annotations;
using EngineLayer.GlycoSearch;

namespace Test
{
    [TestFixture]
    public static class XLTestOGlyco
    {
        [Test]
        public static void OGlycoTest_LoadGlycanBox()
        {
            GlycanBox.GlobalOGlycans = Glycan.LoadGlycan(GlobalVariables.OGlycanLocation).ToArray();
            var GlycanBoxes = GlycanBox.BuildOGlycanBoxes(3);
            Assert.AreEqual(GlycanBoxes.Count(), 454);
        }

        [Test]
        public static void OGlycoTest_GetK()
        {
            List<int> input = new List<int> {1, 2, 3, 4, 5 };

            //Combination test
            var kcombs = Glycan.GetKCombs(input, 3);

            Assert.AreEqual(kcombs.Count(), 10);

            var allcombs = Glycan.GetKCombs(input, 5);

            Assert.AreEqual(allcombs.Count(), 1);

            //Combination test with repetition
            var kcombs_rep = Glycan.GetKCombsWithRept(input, 3);

            Assert.AreEqual(kcombs_rep.Count(), 35);

            var allcombs_rep = Glycan.GetKCombsWithRept(input, 5);

            Assert.AreEqual(allcombs_rep.Count(), 126);

            //Permutation test
            var kperm = Glycan.GetPermutations(input, 3);

            Assert.AreEqual(kperm.Count(), 60);

            var allperm = Glycan.GetPermutations(input, 5).ToList();

            Assert.AreEqual(allperm.Count(), 120);

            //Permutation test with repetition
            var kperm_rep = Glycan.GetPermutationsWithRept(input, 3);

            Assert.AreEqual(kperm_rep.Count(), 125);
        }

        [Test]
        public static void OGlycoTest_GetKPerWithDuplicate()
        {
            List<int> input = new List<int> { 3, 5, 2, 7};
            int[] ids = new int[3]{ 2, 2, 3};
            var perWithDuplicate = GlycoPeptides.GetPermutations(input, ids);
            var allPermutation = Glycan.GetPermutations(input, ids.Length);
            Assert.That(perWithDuplicate.Count() == allPermutation.Count()/2);
        }

        [Test]
        public static void OGlycoTest_StcE()
        {
            Protein p1 = new Protein("MPLFKNTSVSSLYSGTKCR", "");
            var alphaPeptide = p1.Digest(new DigestionParams(protease: "StcE-trypsin"),
                new List<Modification>(), new List<Modification>()).ToArray();
            Assert.That(alphaPeptide.Length == 8);
            Assert.That(alphaPeptide.First().BaseSequence == "MPLFKNTSV");
        }

        [Test]
        public static void OGlycoTest_FragmentIons()
        {
            GlycanBox.GlobalOGlycans = Glycan.LoadGlycan(GlobalVariables.OGlycanLocation).ToArray();
            GlycanBox.GlobalOGlycanModifications = GlycanBox.BuildGlobalOGlycanModifications(GlycanBox.GlobalOGlycans);
            var OGlycanBoxes = GlycanBox.BuildOGlycanBoxes(3).OrderBy(p => p.Mass).ToArray();
            var glycanBox = OGlycanBoxes[8];

            Protein protein = new Protein("PTLFKNVSLYK", "");
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();

            List<int> modPos = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" });

            var peptideWithMod = GlycoPeptides.OGlyGetTheoreticalPeptide(modPos.ToArray(), peptide, OGlycanBoxes[8]);
            Assert.That(peptideWithMod.FullSequence == "PT[O-Glycosylation:H0N1A1G0F0 on X]LFKNVS[O-Glycosylation:H0N1A0G0F0 on X]LYK");

            var fragments_hcd = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.HCD, peptide, peptideWithMod);

            var fragments_ethcd = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.EThcD, peptide, peptideWithMod);

        }
    }
}
