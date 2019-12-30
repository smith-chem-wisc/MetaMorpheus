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

        [Test]
        public static void OGlycoTest_FragmentIonsHash()
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

            var fragments_hcd = peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both);
            var fragmentsMod_hcd = peptideWithMod.Fragment(DissociationType.HCD, FragmentationTerminus.Both);

            var frag_ments_etd = peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both).ToList();
            var fragmentsMod_etd = peptideWithMod.Fragment(DissociationType.ETD, FragmentationTerminus.Both);

            Tuple<int, int[]> keyValuePairs = new Tuple<int, int[]>(8, modPos.ToArray());

            var fragmentsHash_etd = GlycoPeptides.GetFragmentHash(frag_ments_etd, keyValuePairs, OGlycanBoxes, 1000);

            var frag_ments_etd_origin = GlycoPeptides.GetFragmentHash(frag_ments_etd, new Tuple<int, int[]>(0, null), OGlycanBoxes, 1000);

            var overlap = frag_ments_etd_origin.Intersect(fragmentsHash_etd).Count();
        }

        [Test]
        public static void OGlycoTest_AllCombinationsOf()
        {
            List<List<int>> inputs = new List<List<int>>();
            inputs.Add(new List<int> { 1, 2, 3, 4 });
            inputs.Add(new List<int> { 5, 6 });
            inputs.Add(new List<int> { 7, 8 });
            inputs.Add(new List<int> { 9 });
            inputs.Add(new List<int> { 1, 2 });

            var test = ModBox.AllCombinationsOf(inputs.ToArray());
            Assert.That(test.Count == 32);
        }

        [Test]
        public static void OGlycoTest_GetPossibleModSites()
        {
            //Test ModBox.BuildModBoxes
            ModBox.SelectedModifications = new Modification[5];
            ModBox.SelectedModifications[0] = GlobalVariables.AllModsKnownDictionary["Oxidation on M"];
            ModBox.SelectedModifications[1] = GlobalVariables.AllModsKnownDictionary["Hydroxylation on P"];
            ModBox.SelectedModifications[2] = GlobalVariables.AllModsKnownDictionary["Hydroxylation on K"];
            ModBox.SelectedModifications[3] = GlobalVariables.AllModsKnownDictionary["Glucosylgalactosyl on K"];
            ModBox.SelectedModifications[4] = GlobalVariables.AllModsKnownDictionary["Galactosyl on K"];

            var ModBoxes = ModBox.BuildModBoxes(10).Where(p => !p.MotifNeeded.ContainsKey("K") || (p.MotifNeeded.ContainsKey("K") && p.MotifNeeded["K"].Count <= 3)).OrderBy(p => p.Mass).ToArray(); 
            Assert.That(ModBoxes.Length == 860);

            //Test ModBox.GetPossibleModSites
            PeptideWithSetModifications pep = new PeptideWithSetModifications("MGFQGPAGEP[Common Biological:Hydroxylation on P]GPEP[Common Biological:Hydroxylation on P]GQTGPAGAR", GlobalVariables.AllModsKnownDictionary);

            var test = ModBox.GetPossibleModSites(pep, ModBoxes[12]);
            Assert.That(test.Count == 3);

            //Test  ModBox.GetFragmentHash
            PeptideWithSetModifications pep_original = new PeptideWithSetModifications("MGFQGPAGEPGPEPGQTGPAGAR", GlobalVariables.AllModsKnownDictionary);
            var frags = pep_original.Fragment(DissociationType.HCD, FragmentationTerminus.Both);
            var frags_mod = pep.Fragment(DissociationType.HCD, FragmentationTerminus.Both);

            Tuple<int, int>[] tuples = new Tuple<int, int>[2];
            tuples[0] = new Tuple<int, int>(11, 1);
            tuples[1] = new Tuple<int, int>(15, 1);
            var fraghash = ModBox.GetFragmentHash(frags.ToList(), tuples, 1000);
        }

        [Test]
        public static void OGlycoTest_Localization()
        {
            GlycanBox.GlobalOGlycans = Glycan.LoadGlycan(GlobalVariables.OGlycanLocation).ToArray();
            GlycanBox.GlobalOGlycanModifications = GlycanBox.BuildGlobalOGlycanModifications(GlycanBox.GlobalOGlycans);
            var OGlycanBoxes = GlycanBox.BuildOGlycanBoxes(3).OrderBy(p => p.Mass).ToArray();
            var glycanBox = OGlycanBoxes[19];

            Protein protein = new Protein("TTGSLEPSSGASGPQVSSVK", "P16150");
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();
            List<Product> products = peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both).ToList();


            int[] modPos = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" }).OrderBy(p=>p).ToArray();

            var boxes = GlycanBox.BuildOGlycanBoxes(3, glycanBox.GlycanIds).ToArray();

            LocalizationGraph localizationGraph = new LocalizationGraph(modPos.Length, boxes.Length);

            for (int i = 0; i < modPos.Length; i++)
            {
                int maxLength = i + 1;
                int minlength = glycanBox.GlycanIds.Length - (modPos.Length - 1 - i);

                for (int j = 0; j < boxes.Length; j++)
                {
                    if (boxes[j].NumberOfGlycans<= maxLength && boxes[j].NumberOfGlycans >= minlength)
                    {
                        AdjNode adjNode = new AdjNode(i, j, modPos[i], boxes[j]);
                        if (i == 0)
                        {
                            //Get cost
                            adjNode.maxCost = 0;

                        }
                        else
                        {
                            double maxCost = 0;
                            for (int prej = 0; prej <= j; prej++)
                            {
                                if (boxes[j].NumberOfGlycans <= boxes[prej].NumberOfGlycans + 1 &&
                                    (boxes[prej].GlycanIds.Length == 0 || boxes[j].GlycanIds.Any(boxes[prej].GlycanIds.Contains)) &&
                                    localizationGraph.array[i - 1][prej]!=null)
                                {                       
                                    var cost = localizationGraph.array[i-1][prej].maxCost;
                                    if (cost > maxCost)
                                    {
                                        adjNode.Sources.Clear();
                                        adjNode.Costs.Clear();

                                        adjNode.Sources.Add(prej);
                                        adjNode.Costs.Add(cost);
                                        maxCost = cost;
                                    }
                                    else if(cost == maxCost)
                                    {
                                        adjNode.Sources.Add(prej);
                                        adjNode.Costs.Add(cost);
                                    }

                                }
                            }
                            adjNode.maxCost = adjNode.Costs.Max();
                        }

                        localizationGraph.array[i][j] = adjNode;
                    }
                }

            }

        }
    }
}
