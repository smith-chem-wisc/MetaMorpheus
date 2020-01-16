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
    public class XLTestOGlyco
    {
        private static GlycanBox[] OGlycanBoxes { get; set; }

        [OneTimeSetUp]
        public static void Setup()
        {
            GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.GlycanLocations.Where(p => p.Contains("OGlycan.gdb")).First()).ToArray();
            GlycanBox.GlobalOGlycanModifications = GlycanBox.BuildGlobalOGlycanModifications(GlycanBox.GlobalOGlycans);
            OGlycanBoxes = GlycanBox.BuildOGlycanBoxes(3).OrderBy(p => p.Mass).ToArray();
        }

        [Test]
        public static void OGlycoTest_LoadGlycanBox()
        {
            Assert.AreEqual(OGlycanBoxes.Count(), 454);
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
        public static void OxoniumIonAnalysis()
        {
            Assert.That(Glycan.AllOxoniumIons[4] == 13805550);
            Assert.That(Glycan.AllOxoniumIons[5] == 14406607);
            Assert.That(Glycan.AllOxoniumIons[9] == 20408720);
            Assert.That(Glycan.AllOxoniumIons[10] == 27409268);
            Assert.That(Glycan.AllOxoniumIons[12] == 29210324);
            Assert.That(Glycan.AllOxoniumIons[14] == 36614002);

            //Get Scan
            CommonParameters commonParameters = new CommonParameters(dissociationType: DissociationType.EThcD, trimMsMsPeaks: false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_4565.mgf");
            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
            var scan = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).First();

            var productSearchMode = new SinglePpmAroundZeroSearchMode(20);
            var oxoniumIonIntensities = GlycoPeptides.ScanOxoniumIonFilter(scan, productSearchMode, DissociationType.EThcD);

            //Get glycanBox          
            var glycanBox = OGlycanBoxes[19];

            var satifyOxonium = GlycoPeptides.OxoniumIonsAnalysis(oxoniumIonIntensities, glycanBox);
            Assert.That(satifyOxonium);
        
        }

        [Test]
        public static void OGlycoTest_FragmentIons()
        {
            //Get glycanBox
            var glycanBox = OGlycanBoxes[8];

            Protein protein = new Protein("PTLFKNVSLYK", "");
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();

            List<int> modPos = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" });

            var peptideWithMod = GlycoPeptides.OGlyGetTheoreticalPeptide(modPos.ToArray(), peptide, OGlycanBoxes[8]);
            Assert.That(peptideWithMod.FullSequence == "PT[O-Glycosylation:N1A1 on X]LFKNVS[O-Glycosylation:N1 on X]LYK");

            var fragments_hcd = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.HCD, peptide, peptideWithMod);

            var fragments_ethcd = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.EThcD, peptide, peptideWithMod);
        }

        [Test]
        public static void OGlycoTest_FragmentIonsHash()
        {
            //Get glycanBox
            var glycanBox = OGlycanBoxes[8];

            Protein protein = new Protein("PTLFKNVSLYK", "");
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();

            List<int> modPos = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" });

            var peptideWithMod = GlycoPeptides.OGlyGetTheoreticalPeptide(modPos.ToArray(), peptide, glycanBox);
            Assert.That(peptideWithMod.FullSequence == "PT[O-Glycosylation:N1A1 on X]LFKNVS[O-Glycosylation:N1 on X]LYK");

            //The following code prove that the default Fragment method doesn't work for O-glycopeptide due to neutral losses.
            var fragments_hcd = peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both);
            var fragmentsMod_hcd = peptideWithMod.Fragment(DissociationType.HCD, FragmentationTerminus.Both); 
            Assert.That(fragments_hcd.Count() == 20);
            Assert.That(fragmentsMod_hcd.Count() == 50); //The Fragments also contain neutral loss ions. 

            var frag_ments_etd = peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both).ToList();
            var fragmentsMod_etd = peptideWithMod.Fragment(DissociationType.ETD, FragmentationTerminus.Both).ToList();

            //Tuple<int, int[]> keyValuePair represents: <glycanBoxId, glycanModPositions> 
            Tuple<int, int[]> keyValuePairs = new Tuple<int, int[]>(8, modPos.ToArray());

            var fragments_etd_origin = GlycoPeptides.GetFragmentHash(frag_ments_etd, new Tuple<int, int[]>(0, null), OGlycanBoxes, 1000);

            var fragmentsHash_etd = GlycoPeptides.GetFragmentHash(frag_ments_etd, keyValuePairs, OGlycanBoxes, 1000);

            var fragmentsMod_etd_origin = GlycoPeptides.GetFragmentHash(fragmentsMod_etd, new Tuple<int, int[]>(0, null), OGlycanBoxes, 1000);

            var overlap = fragmentsHash_etd.Intersect(fragments_etd_origin).Count();

            Assert.That(overlap == 7);

            var overlap2 = fragmentsHash_etd.Intersect(fragmentsMod_etd_origin).Count();

            Assert.That(overlap2 == 30);
        }

        [Test]
        public static void OGlycoTest_Localization()
        {
            //Get glycanBox
            var glycanBox = OGlycanBoxes[19];

            //Get unmodified peptide, products, allPossible modPos and all boxes.
            Protein protein = new Protein("TTGSLEPSSGASGPQVSSVK", "P16150");
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();
            //var peptide = new PeptideWithSetModifications("TTGSLEPSSGASGPQVSSVK", GlobalVariables.AllModsKnownDictionary);
            List<Product> products = peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both).ToList();

            int[] modPos = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" }).OrderBy(p=>p).ToArray();
            var boxes = GlycanBox.BuildChildOGlycanBoxes(3, glycanBox.ModIds).ToArray();
            Assert.That(boxes.Count() == 6);

            //Test GetLocalFragmentHash, which is used for localiation.
            var testProducts = GlycoPeptides.GetLocalFragmentHash(products, peptide.Length, modPos, 0, glycanBox, boxes[1], 1000);
            var testProducts1 = GlycoPeptides.GetLocalFragmentHash(products, peptide.Length, modPos, 1, glycanBox, boxes[1], 1000);
            Assert.That(testProducts.Count() == 2);
            Assert.That(testProducts1.Count() == 4);

            //Get hashset int
            CommonParameters commonParameters = new CommonParameters(dissociationType: DissociationType.EThcD, trimMsMsPeaks:false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_4565.mgf");
            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
            var scans = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).ToArray();

            int obsPreviousFragmentCeilingMz = 0;
            List<int> binsToSearch = new List<int>();
            foreach (var envelope in scans.First().ExperimentalFragments)
            {
                // assume charge state 1 to calculate mass tolerance
                double experimentalFragmentMass = envelope.monoisotopicMass;

                // get theoretical fragment bins within mass tolerance
                int obsFragmentFloorMass = (int)Math.Floor((commonParameters.ProductMassTolerance.GetMinimumValue(experimentalFragmentMass)) * 1000);
                int obsFragmentCeilingMass = (int)Math.Ceiling((commonParameters.ProductMassTolerance.GetMaximumValue(experimentalFragmentMass)) * 1000);

                // prevents double-counting peaks close in m/z and lower-bound out of range exceptions
                if (obsFragmentFloorMass < obsPreviousFragmentCeilingMz)
                {
                    obsFragmentFloorMass = obsPreviousFragmentCeilingMz;
                }
                obsPreviousFragmentCeilingMz = obsFragmentCeilingMass + 1;

                // search mass bins within a tolerance
                for (int fragmentBin = obsFragmentFloorMass; fragmentBin <= obsFragmentCeilingMass; fragmentBin++)
                {
                    binsToSearch.Add(fragmentBin);
                }
            }
            HashSet<int> allPeaks = new HashSet<int>(binsToSearch);

            //Known peptideWithMod match.
            var peptideWithMod = GlycoPeptides.OGlyGetTheoreticalPeptide(new int[3] { 10, 2, 3}, peptide, glycanBox);
            Assert.That(peptideWithMod.FullSequence == "T[O-Glycosylation:H1N1 on X]T[O-Glycosylation:H1N1 on X]GSLEPSS[O-Glycosylation:N1 on X]GASGPQVSSVK");
            //List<Product> knownProducts = peptideWithMod.Fragment(DissociationType.EThcD, FragmentationTerminus.Both).ToList();
            List<Product> knownProducts = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.EThcD, peptide, peptideWithMod);
            var matchedKnownFragmentIons = MetaMorpheusEngine.MatchFragmentIons(scans.First(), knownProducts, commonParameters);

            //Graph Localization
            LocalizationGraph localizationGraph = new LocalizationGraph(modPos.Length, boxes.Length);

            localizationGraph.LocalizeOGlycan(modPos, glycanBox, boxes, allPeaks, products, peptide.Length);

            var allPaths = LocalizationGraph.GetAllPaths(localizationGraph.array, boxes);

            var knowPath = new int[8] {2, 4, 4, 4, 5, 5, 5, 5 };
            Assert.That(Enumerable.SequenceEqual(knowPath, allPaths[0]));

            var local = LocalizationGraph.GetLocalizedPath(localizationGraph.array, modPos, boxes, allPaths.First());

            var knowLocal = new Tuple<int, int>[3] { new Tuple<int, int>(2, 1), new Tuple<int, int>(3, 1), new Tuple<int, int>(10, 0) };
            Assert.That(Enumerable.SequenceEqual(local, knowLocal));
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

            //Test  ModBox.GetFragmentHash
            PeptideWithSetModifications pep_original = new PeptideWithSetModifications("MGFQGPAGEPGPEPGQTGPAGAR", GlobalVariables.AllModsKnownDictionary);
            var frags = pep_original.Fragment(DissociationType.HCD, FragmentationTerminus.Both);
            var frags_mod = pep.Fragment(DissociationType.HCD, FragmentationTerminus.Both);

            Tuple<int, int>[] tuples = new Tuple<int, int>[2];
            tuples[0] = new Tuple<int, int>(11, 1);
            tuples[1] = new Tuple<int, int>(15, 1);

        }

        [Test]
        public static void OGlycoTest_GetLeft()
        {
            int[] array1 = new int[6] { 0, 0, 0, 1, 1, 2 };
            int[] array2 = new int[3] { 0, 0, 1 };
            var left = LocalizationGraph.GetLeft(array1, array2);

            var knowLeft = new int[3] { 0, 1, 2 };
            Assert.That(Enumerable.SequenceEqual(left, knowLeft));
        }

        [Test]
        public static void OGlycoTest_LocalizationMod()
        {
            //Get modBox
            ModBox.SelectedModifications = new Modification[4];
            ModBox.SelectedModifications[0] = GlobalVariables.AllModsKnownDictionary["Hydroxylation on P"];
            ModBox.SelectedModifications[1] = GlobalVariables.AllModsKnownDictionary["Hydroxylation on K"];
            ModBox.SelectedModifications[2] = GlobalVariables.AllModsKnownDictionary["Glucosylgalactosyl on K"];
            ModBox.SelectedModifications[3] = GlobalVariables.AllModsKnownDictionary["Galactosyl on K"];

            var ModBoxes = ModBox.BuildModBoxes(10).Where(p => !p.MotifNeeded.ContainsKey("K") || (p.MotifNeeded.ContainsKey("K") && p.MotifNeeded["K"].Count <= 3)).OrderBy(p => p.Mass).ToArray();
            var modBox = ModBoxes[10];

            //Get unmodified peptide, products, allPossible modPos and all boxes.
            Protein protein = new Protein("GSDGSVGPVGPAGPIGSAGPPGFPGAPGPKGEIGAVGNAGPAGPAGPR", "P08123");
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).ElementAt(2);
            //PeptideWithSetModifications peptide = new PeptideWithSetModifications("GSDGSVGPVGPAGPIGSAGPPGFPGAPGPKGEIGAVGNAGPAGPAGPR", GlobalVariables.AllModsKnownDictionary);
            var products = peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both).ToList();
            int[] modPos = ModBox.GetAllPossibleModSites(peptide, modBox);
            var boxes = ModBox.BuildChildModBoxes(modBox.NumberOfMods, modBox.ModIds).ToArray();
            Assert.That(boxes.Count() == 8);
            //Test boxSatisfyBox
            var boxSatisfyBox = LocalizationGraph.BoxSatisfyBox(boxes);

            //Test GetAllPossibleModSites
            var testModPos = ModBox.GetAllPossibleModSites(peptide, ModBoxes[11]);
            Assert.That(testModPos == null);

            //Test GetLocalFragmentHash, which is used for localiation.
            var testProducts = ModBox.GetLocalFragmentHash(products, peptide.Length, modPos, 0, modBox, boxes[1], 1000);
            var testProducts1 = ModBox.GetLocalFragmentHash(products, peptide.Length, modPos, 2, modBox, boxes[1], 1000);
            Assert.That(testProducts.Count() == 6);
            Assert.That(testProducts1.Count() == 12);

            //Test BoxSatisfyModPos
            var pos1 = LocalizationGraph.BoxSatisfyModPos(modBox, boxes[5], 30, peptide);
            Assert.That(pos1);
            var pos2 = LocalizationGraph.BoxSatisfyModPos(modBox, boxes[5], 31, peptide);
            Assert.That(!pos2);

            //Get hashset int
            CommonParameters commonParameters = new CommonParameters(dissociationType: DissociationType.HCD, trimMsMsPeaks: false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\id_08-24-19_AC-P_patient146_0p12mg_0p9uL-calib_20086.mgf");
            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
            var scans = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).ToArray();

            int obsPreviousFragmentCeilingMz = 0;
            List<int> binsToSearch = new List<int>();
            foreach (var envelope in scans.First().ExperimentalFragments)
            {
                // assume charge state 1 to calculate mass tolerance
                double experimentalFragmentMass = envelope.monoisotopicMass;

                // get theoretical fragment bins within mass tolerance
                int obsFragmentFloorMass = (int)Math.Floor((commonParameters.ProductMassTolerance.GetMinimumValue(experimentalFragmentMass)) * 1000);
                int obsFragmentCeilingMass = (int)Math.Ceiling((commonParameters.ProductMassTolerance.GetMaximumValue(experimentalFragmentMass)) * 1000);

                // prevents double-counting peaks close in m/z and lower-bound out of range exceptions
                if (obsFragmentFloorMass < obsPreviousFragmentCeilingMz)
                {
                    obsFragmentFloorMass = obsPreviousFragmentCeilingMz;
                }
                obsPreviousFragmentCeilingMz = obsFragmentCeilingMass + 1;

                // search mass bins within a tolerance
                for (int fragmentBin = obsFragmentFloorMass; fragmentBin <= obsFragmentCeilingMass; fragmentBin++)
                {
                    binsToSearch.Add(fragmentBin);
                }
            }
            HashSet<int> allPeaks = new HashSet<int>(binsToSearch);

            //Known peptideWithMod match.
            var peptideWithMod = ModBox.GetTheoreticalPeptide(new int[4] {22, 25, 28, 31}, peptide, modBox);
            Assert.That(peptideWithMod.FullSequence == "GSDGSVGPVGPAGPIGSAGPP[Common Biological:Hydroxylation on P]GFP[Common Biological:Hydroxylation on P]GAP[Common Biological:Hydroxylation on P]GPK[Common Biological:Hydroxylation on K]GEIGAVGNAGPAGPAGPR");
            List<Product> knownProducts = peptideWithMod.Fragment(DissociationType.HCD, FragmentationTerminus.Both).ToList();
            var matchedKnownFragmentIons = MetaMorpheusEngine.MatchFragmentIons(scans.First(), knownProducts, commonParameters);
            Assert.That(matchedKnownFragmentIons.Count == 42);

            //Graph Localization
            LocalizationGraph localizationGraph = new LocalizationGraph(modPos.Length, boxes.Length);
            localizationGraph.LocalizeMod(modPos, modBox, boxes, allPeaks, products, peptide);

            var allPaths = LocalizationGraph.GetAllPaths(localizationGraph.array, boxes);

            var knowPath = new int[12] { 0, 0, 0, 0, 1, 3, 5, 5, 7, 7, 7, 7 };
            Assert.That(Enumerable.SequenceEqual(knowPath, allPaths[0]));

            var local = LocalizationGraph.GetLocalizedPath(localizationGraph.array, modPos, boxes, allPaths.First());

            var knowLocal = new Tuple<int, int>[4] {new Tuple<int, int>(22, 0), new Tuple<int, int>(25, 0) , new Tuple<int, int>(28, 0) , new Tuple<int, int>(31, 1) };
            Assert.That(Enumerable.SequenceEqual(local, knowLocal));
        }

        [Test]
        public static void OGlcoTest_GetAllPaths()
        {
            int[] modPos = new int[3] { 2, 4, 6 };
            var glycanBox = OGlycanBoxes[19];
            var boxes = GlycanBox.BuildChildOGlycanBoxes(3, glycanBox.ModIds).ToArray();
            LocalizationGraph localizationGraph = new LocalizationGraph(modPos.Length, boxes.Length);

            for (int i = 0; i < modPos.Length; i++)
            {
                for (int j = 0; j < boxes.Length; j++)
                {
                    localizationGraph.array[i][j] = new AdjNode(i, j, modPos[i], boxes[j]);
                    localizationGraph.array[i][j].Sources = new List<int> { j };
                    localizationGraph.array[i][j].Costs = new List<double> { 1 };
                    localizationGraph.array[i][j].maxCost = 1;
                }
            }
            localizationGraph.array[2][5].Sources = new List<int> {  4, 5 };
            localizationGraph.array[2][5].Costs = new List<double> { 1, 1 };

            localizationGraph.array[1][4].Sources = new List<int> { 2, 4 };
            localizationGraph.array[1][4].Costs = new List<double> { 1, 1 };

            var allPaths = LocalizationGraph.GetAllPaths(localizationGraph.array, boxes);

            Assert.That(allPaths.Count == 3);
            Assert.That(allPaths[0] == new int[3] { 2, 4, 5});
            Assert.That(allPaths[1] == new int[3] { 4, 4, 5 });
            Assert.That(allPaths[2] == new int[3] { 5, 5, 5 });
        }
    }
}
