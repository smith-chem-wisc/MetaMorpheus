using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;
using Nett;
using EngineLayer.GlycoSearch;
using FlashLFQ;
using NUnit.Framework.Internal;
using SpectralAveraging;
using Chemistry;
using MzLibUtil;
using Readers;
using System.Text;
using Omics.Modifications;
using ThermoFisher.CommonCore.BackgroundSubtraction;
using Easy.Common.Extensions;
using iText.IO.Font.Otf;
using static Nett.TomlObjectFactory;
using Omics.SpectrumMatch;
using TopDownProteomics;

namespace Test
{
    [TestFixture]
    public class TestOGlyco
    {
        private static GlycanBox[] OGlycanBoxes { get; set; }


        [OneTimeSetUp]
        public static void Setup()
        {
            GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.OGlycanLocations.Where(p => p.Contains("OGlycan.gdb")).First(), true, true).ToArray();
            GlycanBox.GlobalOGlycanModifications = GlycanBox.BuildGlobalOGlycanModifications(GlycanBox.GlobalOGlycans);
            OGlycanBoxes = GlycanBox.BuildOGlycanBoxes(3).OrderBy(p => p.Mass).ToArray();
        }

        [Test]
        public static void OGlycoTest_LoadGlycanBox()
        {
            Assert.AreEqual(OGlycanBoxes.Count(), 454);
        }

        [Test]
        public static void OGlycanTest_GetGlycanBox_Decoy()
        {
            GlycanBox[] OGlycanBoxes = GlycanBox.BuildOGlycanBoxes(3).ToArray();
            Assert.That(OGlycanBoxes.All(p => p.TargetDecoy = true));

            GlycanBox[] OGlycanBoxes_withDecoys = GlycanBox.BuildOGlycanBoxes(3, true).ToArray();
            var group_target = OGlycanBoxes_withDecoys.GroupBy(p => p.TargetDecoy == true);
            var group_decoy = OGlycanBoxes_withDecoys.GroupBy(p => p.TargetDecoy == false);
            Assert.That(group_target.Count() == group_decoy.Count());

        }

        [Test]
        public static void GlycoTest_WritingSummary()  // In order to test writing function on different search type ex. O-Search, N-search, N-O search, make sure we have the corresponding search Rseult file.
        {
            string outputFolder_NSearch = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder_NSearch);

            var glycoSearchTask_NSearch = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSearchTaskconfigNGlycoTest_Run.toml"), MetaMorpheusTask.tomlConfig);

            DbForTask db = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P16150.fasta"), false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask_NSearch) }, new List<string> { spectraFile }, new List<DbForTask> { db }, outputFolder_NSearch).Run();

            Directory.Delete(outputFolder_NSearch, true);


            string outputFolder_NOSearch = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder_NOSearch);

            var glycoSearchTask_NOSearch = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSearchTaskconfigN_OGlycoTest_Run.toml"), MetaMorpheusTask.tomlConfig);

            string spectraFile_NOSearch = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask_NOSearch) }, new List<string> { spectraFile_NOSearch }, new List<DbForTask> { db }, outputFolder_NOSearch).Run();

            Directory.Delete(outputFolder_NOSearch, true);
        }


        [Test]
        public static void OGlycanTest_IsobaricCase()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\Glyco_Isobaric_testing.toml"), MetaMorpheusTask.tomlConfig);
            glycoSearchTask._glycoSearchParameters.DoParsimony = false;
            glycoSearchTask._glycoSearchParameters.DoQuantification = true;
            glycoSearchTask._glycoSearchParameters.OxoniumIonFilt = false; //turn off the diagnostic filter

            DbForTask targetDbForTask = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoProteinFASTA_7proteins.fasta"), false);
            DbForTask contaminDbForTask = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P13987_contaminant.fasta"), true);

            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoPepMix_snip.mzML");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) }, new List<string> { spectraFile }, new List<DbForTask> { targetDbForTask, contaminDbForTask }, outputFolder).Run();
            
            string oGlycoPath = Path.Combine(outputFolder, "Task", "oglyco.psmtsv");            
            var glycanLevel_filterOFF = PsmTsvReader.ReadTsv(oGlycoPath, out var errors) //load the PSMs data from the "csv file" and bulid the objects
                .Where(p => p.Ms2ScanNumber == 161 && p.BaseSeq == "HTSVQTTSSGSGPFTDVR").ToList()[0].GlycanLocalizationLevel;

            
            Assert.That(glycanLevel_filterOFF != EngineLayer.GlycoSearch.LocalizationLevel.Level1 && glycanLevel_filterOFF != EngineLayer.GlycoSearch.LocalizationLevel.Level1b);
            Directory.Delete(outputFolder, true);

            glycoSearchTask._glycoSearchParameters.OxoniumIonFilt = true;
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) }, new List<string> { spectraFile }, new List<DbForTask> { targetDbForTask, contaminDbForTask }, outputFolder).Run();
            var glycanLevel_filterON = PsmTsvReader.ReadTsv(oGlycoPath, out var error) //load the PSMs data from the "csv file" and bulid the objects
                .Where(p => p.Ms2ScanNumber == 161 && p.BaseSeq == "HTSVQTTSSGSGPFTDVR").ToList()[0].GlycanLocalizationLevel;


            Assert.That(glycanLevel_filterON == EngineLayer.GlycoSearch.LocalizationLevel.Level1);
            Directory.Delete(outputFolder, true);

        }



        [Test]
        public static void GlycoSpectralHeader()
        {
            string header = GlycoSpectralMatch.GetTabSepHeaderSingle().Trim();
            string expectedHeader = "File Name\tScan Number\tRetention Time\tPrecursor Scan Number\tPrecursor MZ\tPrecursor Charge\tPrecursor Mass\tProtein Accession\tOrganism\tProtein Name\tStart and End Residues In Protein\tBase Sequence\tFlankingResidues\tFull Sequence\tNumber of Mods\tPeptide Monoisotopic Mass\tScore\tRank\tMatched Ion Series\tMatched Ion Mass-To-Charge Ratios\tMatched Ion Mass Diff (Da)\tMatched Ion Mass Diff (Ppm)\tMatched Ion Intensities\tMatched Ion Counts\tDecoy/Contaminant/Target\tQValue\tPEP\tPEP_QValue\t";
            Assert.That(header.Equals(expectedHeader.Trim()));
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
        public static void OGlycoTest_OGlycanChildIons()
        {
            // Reload the glycan database to test the child ions.
            GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.OGlycanLocations.Where(p => p.Contains("OGlycan.gdb")).First(), true, true).ToArray();

            var glycan = GlycanBox.GlobalOGlycans[5]; // we use the glycan (N(H)(N(H)))

            Assert.That(glycan.Ions.Count == 5);

            var kind = glycan.Kind;

            var glycanIons = GlycanDatabase.OGlycanCompositionCombinationChildIons(kind);

            Assert.That(glycanIons.Count() == 6);

            var coreIons = GlycanDatabase.OGlycanCompositionFragments(kind);
            Assert.That(coreIons.Count() == 6);

            //The following code is to test the glycan with complex structure, only to pass the converage.

            var testKind = GlycanDatabase.String2Kind("HexNAc(2)Hex(4)Fuc(2)NeuAc(1)Xylose(1)");
            
            var testGlycanIons = GlycanDatabase.OGlycanCompositionFragments(testKind);


            var testKind_smallGlycan = GlycanDatabase.String2Kind("HexNAc(1)");

            var testGlycanIons_smallGlycan = GlycanDatabase.OGlycanCompositionFragments(testKind_smallGlycan);


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
        public static void OGlycanTest_Kdn()
        {
            var kind = new byte[] { 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2 };
            double mass = Glycan.GetMass(kind) / 1E5;
            string name = Glycan.GetKindString(kind);
            Assert.AreEqual(name, "H1N1K2");
            Assert.AreEqual(mass, 865.27013);

            string kdnGlycan = "HexNAc(2)Hex(2)Kdn(1)";
            string kdnGlycan2 = "N(H)H(N)K";
            var expectedKind = new byte[] { 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1 };
            Assert.AreEqual(GlycanDatabase.String2Kind(kdnGlycan), expectedKind);
            Assert.AreEqual(Glycan.GetKind(kdnGlycan2), expectedKind);
        }

        [Test]
        public static void GlycoTest_MotifExist()
        {
            string baseSeq = "AFGQFFSPGEVIYNKTDRAG";
            var exist = GlycoSpectralMatch.MotifExist(baseSeq, new string[] { "Nxt", "Nxs" });
            Assert.That(exist);
        }

        [Test] // In this test, there are 272 and 294 oxonium ions in the scan, but the glycanBox doesn't contain these HexNAc ions.
        public static void DiagonsticFilter()
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
            var oxoniumIonIntensities = GlycoPeptides.ScanOxoniumIonFilter(scan, productSearchMode);

            //Get glycanBox          
            var glycanBox = OGlycanBoxes[19];

            var satifyOxonium = GlycoPeptides.DiagonsticFilter(oxoniumIonIntensities, glycanBox);
            Assert.That(!satifyOxonium);
        
        }

        [Test]
        public static void OGlycoTest_FragmentIons()
        {
            //This test is to test ETD on proline. As proline didn't have proline related product ions.

            //Get glycanBox
            var glycanBox = OGlycanBoxes[8];

            Protein protein = new Protein("PTLFKNVSLYK", "");
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();

            List<int> modPos = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" });

            var peptideWithMod = GlycoPeptides.OGlyGetTheoreticalPeptide(modPos.ToArray(), peptide, OGlycanBoxes[8]);
            Assert.That(peptideWithMod.FullSequence == "PT[O-Glycosylation:N1A1 on X]LFKNVS[O-Glycosylation:N1 on X]LYK");

            var fragments_hcd = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.HCD, new List<ProductType>(), peptide, peptideWithMod);

            var fragments_ethcd = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.EThcD, new List<ProductType>(), peptide, peptideWithMod);

            var customIons = new List<ProductType> { ProductType.c, ProductType.zDot };
            var fragments_custom = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.Custom, customIons, peptide, peptideWithMod);

            Assert.That(fragments_hcd.Count == fragments_custom.Count && fragments_hcd.Count == fragments_ethcd.Count - 20);
        }

        [Test]
        public static void OGlycoTest_FragmentIons2()
        {
            //Get glycanBox
            var glycanBox = OGlycanBoxes[24];
            Protein protein = new Protein("TVYLGASK", "");
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();

            List<int> modPos = new List<int> { 2, 8 };

            var peptideWithMod = GlycoPeptides.OGlyGetTheoreticalPeptide(modPos.ToArray(), peptide, glycanBox);
            Assert.That(peptideWithMod.FullSequence == "T[O-Glycosylation:H1N1 on X]VYLGAS[O-Glycosylation:H1N1A1 on X]K");

            var fragments_etd = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.ETD, new List<ProductType>(), peptide, peptideWithMod);


            Assert.That(fragments_etd.Count == 22);
            Assert.That(fragments_etd.Last().Annotation == "zDot8");
            Assert.That(fragments_etd.Last().NeutralMass > 1824);

            var fragments_hcd = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.HCD, new List<ProductType>(), peptide, peptideWithMod);
            Assert.That(fragments_hcd.Where(p=>p.ProductType == ProductType.M).Count() == 5);

            var fragments_ethcd = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.EThcD, new List<ProductType>(), peptide, peptideWithMod);
            Assert.That(fragments_ethcd.Where(p => p.ProductType == ProductType.M).Count() == 5);
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
            var fragments_hcd = new List<Product>();
                peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, fragments_hcd);
            var fragmentsMod_hcd = new List<Product>();
                peptideWithMod.Fragment(DissociationType.HCD, FragmentationTerminus.Both, fragmentsMod_hcd); 
            Assert.That(fragments_hcd.Count() == 20);
            Assert.That(fragmentsMod_hcd.Count() == 61); //The Fragments also contain neutral loss ions. 

            var frag_ments_etd = new List<Product>();
                peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, frag_ments_etd);
            var fragmentsMod_etd = new List<Product>();
                peptideWithMod.Fragment(DissociationType.ETD, FragmentationTerminus.Both, fragmentsMod_etd);

            //Tuple<int, int[]> keyValuePair represents: <glycanBoxId, glycanModPositions> 
            Tuple<int, int[]> keyValuePairs = new Tuple<int, int[]>(8, modPos.ToArray());

            var fragments_etd_origin = GlycoPeptides.GetFragmentHash(frag_ments_etd, new Tuple<int, int[]>(0, null), OGlycanBoxes, 1000);

            var fragmentsHash_etd = GlycoPeptides.GetFragmentHash(frag_ments_etd, keyValuePairs, OGlycanBoxes, 1000);

            var fragmentsMod_etd_origin = GlycoPeptides.GetFragmentHash(fragmentsMod_etd, new Tuple<int, int[]>(0, null), OGlycanBoxes, 1000);

            var overlap = fragmentsHash_etd.Intersect(fragments_etd_origin).Count();

            Assert.That(overlap == 14);

            var overlap2 = fragmentsHash_etd.Intersect(fragmentsMod_etd_origin).Count();

            //ETD didn't change y ions.
            Assert.That(overlap2 == 23);
        }

        [Test]
        public static void OGlycoTest_Localization()
        {
            //Get glycanBox
            var glycanBox = OGlycanBoxes[19];

            //Get unmodified peptide, products, allPossible modPos and all boxes.
            Protein protein = new Protein("TTGSLEPSSGASGPQVSSVK", "P16150");
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();
            List<Product> products = new List<Product>();
            peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);

            int[] modPos = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" }).OrderBy(v=>v).ToArray();
            var boxes = GlycanBox.BuildChildOGlycanBoxes(3, glycanBox.ModIds).ToArray();
            Assert.That(boxes.Count() == 6);

            //Get Unlocal Fragment
            var unlocalCost = GlycoPeptides.GetUnlocalFragment(products, modPos, glycanBox);
            Assert.That(unlocalCost.Count == 4); //Basicly, the unlocal are c/z ions that don't localize glycosylation. 

            //Get scan
            CommonParameters commonParameters = new CommonParameters(dissociationType: DissociationType.EThcD, trimMsMsPeaks:false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_4565.mgf");
            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
            var scans = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).ToArray();

            //Known peptideWithMod match.
            var peptideWithMod = GlycoPeptides.OGlyGetTheoreticalPeptide(new int[3] { 10, 2, 3}, peptide, glycanBox);
            Assert.That(peptideWithMod.FullSequence == "T[O-Glycosylation:H1N1 on X]T[O-Glycosylation:H1N1 on X]GSLEPSS[O-Glycosylation:N1 on X]GASGPQVSSVK");
            List<Product> knownProducts = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.EThcD, new List<ProductType>(), peptide, peptideWithMod);
            var matchedKnownFragmentIons = MetaMorpheusEngine.MatchFragmentIons(scans.First(), knownProducts, commonParameters);

            //Graph Localization
            LocalizationGraph localizationGraph = new LocalizationGraph(modPos, glycanBox, boxes, -1);
            LocalizationGraph.LocalizeOGlycan(localizationGraph, scans.First(), commonParameters.ProductMassTolerance, products);
            var allPaths = LocalizationGraph.GetAllHighestScorePaths(localizationGraph.array, localizationGraph.ChildModBoxes);
            var knowPath = new int[8] {2, 4, 4, 4, 5, 5, 5, 5 };
            Assert.That(Enumerable.SequenceEqual(knowPath, allPaths[0]));

            //Get localized Route
            var local = LocalizationGraph.GetLocalizedPath(localizationGraph, allPaths.First());
            Assert.That(Enumerable.SequenceEqual(local.Mods.Select(v=>v.Item1), new List<int>{ 2, 3, 10}));
            Assert.That(Enumerable.SequenceEqual(local.Mods.Select(v => v.Item2), new List<int> { 1, 1, 0 }));


            //Get all paths, calculate PScore and calculate position probability. 
            var p = scans.First().TheScan.MassSpectrum.Size * commonParameters.ProductMassTolerance.GetRange(1000).Width / scans.First().TheScan.MassSpectrum.Range.Width;
            var n = knownProducts.Where(v=>v.ProductType == ProductType.c || v.ProductType == ProductType.zDot).Count();
            var allPathWithWeights = LocalizationGraph.GetAllPaths_CalP(localizationGraph, p, n);
            Assert.That(allPathWithWeights.Count == 168);

            //Calculate Site Specific Localization Probability
            var y = LocalizationGraph.CalSiteSpecificLocalizationProbability(allPathWithWeights, localizationGraph.ModPos);
            Assert.That(y.Count == 8);
            Assert.That(y.First().Value[1].Item2 > 0.99);

        }

        [Test]
        public static void OGlycoTest_Localization2()
        {
            //There may have a bug that MM cannot identify Peptide modified with (HexNAc), This is to test and find the bug.
            //Get glycanBox
            var glycanBox = OGlycanBoxes[0];

            //Get unmodified peptide, products, allPossible modPos and all boxes.
            Protein protein = new Protein("AATVGSLAGQPLQER", "P16150");
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();
            List<Product> products = new List<Product>();
            peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);

            int[] modPos = GlycoSpectralMatch.GetPossibleModSites(peptide, new string[] { "S", "T" }).OrderBy(p => p).ToArray();
            var boxes = GlycanBox.BuildChildOGlycanBoxes(glycanBox.NumberOfMods, glycanBox.ModIds).ToArray();

            //Load scan.
            CommonParameters commonParameters = new CommonParameters(dissociationType: DissociationType.ETD, trimMsMsPeaks: false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\181217_Fusion_(LC2)_NewObj_Serum_deSA_Jacalin_HRM_4h_ETD_HCD_DDA_mz(400_1200)_21707.mgf");
            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
            var scans = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).ToArray();

            //Known peptideWithMod match.
            var peptideWithMod = GlycoPeptides.OGlyGetTheoreticalPeptide(new int[1] {4}, peptide, glycanBox);
            Assert.That(peptideWithMod.FullSequence == "AAT[O-Glycosylation:N1 on X]VGSLAGQPLQER");
            //List<Product> knownProducts = peptideWithMod.Fragment(DissociationType.EThcD, FragmentationTerminus.Both).ToList();
            List<Product> knownProducts = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.ETD, new List<ProductType>(), peptide, peptideWithMod);
            var matchedKnownFragmentIons = MetaMorpheusEngine.MatchFragmentIons(scans.First(), knownProducts, commonParameters);

            //Get hashset int
            int obsPreviousFragmentCeilingMz = 0;
            List<int> binsToSearch = new List<int>();
            foreach (var envelope in scans.First().ExperimentalFragments)
            {
                // assume charge state 1 to calculate mass tolerance
                double experimentalFragmentMass = envelope.MonoisotopicMass;

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


            //Graph Localization
            LocalizationGraph localizationGraph = new LocalizationGraph(modPos, glycanBox, boxes, -1);

            LocalizationGraph.LocalizeOGlycan(localizationGraph, scans.First(), commonParameters.ProductMassTolerance, products);

            var allPaths = LocalizationGraph.GetAllHighestScorePaths(localizationGraph.array, localizationGraph.ChildModBoxes);

            var knowPath = new int[2] { 1, 1 };
            Assert.That(Enumerable.SequenceEqual(knowPath, allPaths[0]));

            var local = LocalizationGraph.GetLocalizedPath(localizationGraph, allPaths.First());

            Assert.That(Enumerable.SequenceEqual(local.Mods.Select(p=>p.Item1), new List<int> { 4}));
            Assert.That(Enumerable.SequenceEqual(local.Mods.Select(p => p.Item2), new List<int> { 0 }));
        }

        [Test]
        public static void OGlycoTest_Run()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSearchTaskconfigOGlycoTest_Run.toml"), MetaMorpheusTask.tomlConfig);

            DbForTask db = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P16150.fasta"), false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) }, new List<string> { spectraFile }, new List<DbForTask> { db }, outputFolder).Run();

            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void OGlycoTest_Run2()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSearchTaskconfigOGlycoTest_Run2.toml"), MetaMorpheusTask.tomlConfig);

            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P02649.fasta"), false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\181217_Fusion_(LC2)_NewObj_Serum_deSA_Jacalin_HRM_4h_ETD_HCD_DDA_mz(400_1200)_21707.mgf");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) }, new List<string> { spectraFile }, new List<DbForTask> { db }, outputFolder).Run();
            Directory.Delete(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"), true);
        }

        [Test]
        public static void OGlycoTest_Run3()
        {
            //Test search EThcD scan with setting DissociationType ETD.The Oxonium ion always matches, but the filter is user controlled.
            //The scan is modified to contain no 274 and 292 oxonium ions.
            var task = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSearchTaskconfig_ETD_Run3.toml"), MetaMorpheusTask.tomlConfig);

            Directory.CreateDirectory(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"));
            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\FiveMucinFasta.fasta"), false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\22-12-14_EclipseOglyco_EThcD_150ms_calRxn_17360.mgf");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) }, new List<string> { spectraFile }, new List<DbForTask> { db }, Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData")).Run();
            var resultsPath = File.ReadAllLines(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData\Task\oglyco.psmtsv"));
            Assert.That(resultsPath.Length > 1);
            Directory.Delete(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"), true);

            Directory.CreateDirectory(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"));
            var task2 = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSearchTaskconfig_ETD_Run3.toml"), MetaMorpheusTask.tomlConfig);
            task2._glycoSearchParameters.OxoniumIonFilt = true; //turn on the diagnostic filter
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task2) }, new List<string> { spectraFile }, new List<DbForTask> { db }, Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData")).Run();
            var resultsExist = File.Exists(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData\Task\oglyco.psmtsv"));
            Assert.That(!resultsExist);
            Directory.Delete(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"), true);
        }

        //make sure that hyphens in protein names don't produce a crash during protein inference from glycopeptides
        [Test]
        public static void OGlycoTest_Run4()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSearchTaskconfigOGlycoTest_Run.toml"), MetaMorpheusTask.tomlConfig);

            DbForTask db = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P16150withHyphenInName.fasta"), false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) }, new List<string> { spectraFile }, new List<DbForTask> { db }, outputFolder).Run();

            var folders = Directory.GetDirectories(outputFolder).Select(b => Path.GetFileName(b)).ToList();
            var folderContents = Directory.GetFiles(Path.Combine(outputFolder, folders[0])).ToList();
            Assert.That(folderContents[6].Contains("_AllProteinGroups.tsv"));

            Directory.Delete(outputFolder, true);
        }


        [Test]
        public static void OGlycoTest_Run5()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSnip.toml"), MetaMorpheusTask.tomlConfig);
            glycoSearchTask._glycoSearchParameters.WriteContaminants = false;

            DbForTask targetDbForTask = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoProteinFASTA_7proteins.fasta"), false);
            DbForTask contaminDbForTask = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P13987_contaminant.fasta"), true);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoPepMix_snip.mzML");
            new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) },
                new List<string> { spectraFile },
                new List<DbForTask> { targetDbForTask, contaminDbForTask }, 
                outputFolder).Run();


            // TODO: Test output, make sure the values on the results.txt really reflect the number counted in the csv files
            // Parse values from results.txt
            string resultsTextPath = Directory.GetFiles(outputFolder, "allResults.txt", SearchOption.TopDirectoryOnly)
                .FirstOrDefault(); // Try to find the file name "allResults.txt" in the output folder
            if (resultsTextPath is null)
                Assert.Fail("Results file not found.");
            string[] allResultTxtLines = File.ReadAllLines(resultsTextPath); //read all lines from the file
            Assert.That(allResultTxtLines.Length > 0); // make sure there are lines in the file

            //For PSMs
            var allPsmLine = allResultTxtLines.First(p => p.Contains("target PSMs within"));
            int psmCount = int.Parse(allPsmLine.Split(':').Last().Trim());

            //For ProteinGroups
            var proteinGroupLine = allResultTxtLines.First(p => p.Contains("protein groups within"));
            int proteinGroupCount = int.Parse(proteinGroupLine.Split(':').Last().Trim());

            //For GlycoPSMs
            var glycoPsmLine = allResultTxtLines.First(p => p.Contains("O-Glyco PSMs within"));
            int glycoPsmCount = int.Parse(glycoPsmLine.Split(':').Last().Trim()); // read the number of glyco PSMs from the results file

            //For Level1GlycoPSMs
            var level1PsmLine = allResultTxtLines.First(p => p.Contains("Level 1 O-Glyco PSMs within"));
            int level1Psmcount = int.Parse(level1PsmLine.Split(':').Last().Trim()); // read the number of Level1-PSMs from the results file

            // Parse counted number from csv files

            //For PSMs
            var allPsmPath = Path.Combine(outputFolder, "Task", "AllPSMs.psmtsv");
            List<PsmFromTsv> onePercentPsms1 = PsmTsvReader.ReadTsv(allPsmPath, out var errors2)
            .Where(p => p.QValue <= 0.01).ToList();
            Assert.That(errors2.Count == 0);// if we cannot find the file, we will get an error message
            int readInPsmsCount = onePercentPsms1.Count;

            //For ProteinGroups
            var allProteinGroupsPath = Path.Combine(outputFolder, "Task", "_AllProteinGroups.tsv");
            string[] proteinGroupHeaders = File.ReadAllLines(allProteinGroupsPath).First().Split("\t");
            int readInProteinCount = File.ReadAllLines(Path.Combine(outputFolder, "Task", "_AllProteinGroups.tsv")).Skip(1)
                .Select(line => line.Split('\t'))
                .Count(p => double.TryParse(p[Array.IndexOf(proteinGroupHeaders, "Protein QValue")], out double qVaule) && qVaule < 0.01);

            //For GlycoPSMs
            string oGlycoPath = Path.Combine(outputFolder, "Task", "oglyco.psmtsv");
            List<PsmFromTsv> onePercentoGlycoPsms = PsmTsvReader.ReadTsv(oGlycoPath, out var errors) //load the PSMs data from the "csv file" and bulid the objects
                .Where(p => p.QValue <= 0.01).ToList(); // the filtering (Q<0.01)
            int readInGlycoPsmCount = onePercentoGlycoPsms.Count; // the gPSMs number with Fdr<0.01
            Assert.That(errors.Count == 0);// if we cannot find the file, we will get an error message

            //For Level1GlycoPSMs
            int readInLevel1GlycoPsmCount = onePercentoGlycoPsms.Count(p => p.GlycanLocalizationLevel == EngineLayer.GlycoSearch.LocalizationLevel.Level1); //the level1 gPSMs number

            //Compare the numbers
            Assert.That(psmCount, Is.EqualTo(readInPsmsCount));
            Assert.That(proteinGroupCount, Is.EqualTo(readInProteinCount));
            Assert.That(glycoPsmCount, Is.EqualTo(readInGlycoPsmCount));
            Assert.That(level1Psmcount, Is.EqualTo(readInLevel1GlycoPsmCount));


            Directory.Delete(outputFolder, true);
        }

        [Test] 
        public static void OGlycoTest_Run5_WriteContaminants()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSnip.toml"), MetaMorpheusTask.tomlConfig);
            glycoSearchTask._glycoSearchParameters.WriteContaminants = true; // write contaminants to the output folder

            DbForTask targetDbForTask = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoProteinFASTA_7proteins.fasta"), false);
            DbForTask contaminDbForTask = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P13987_contaminant.fasta"), true);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoPepMix_snip.mzML");
            new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) },
                new List<string> { spectraFile },
                new List<DbForTask> { targetDbForTask, contaminDbForTask },
                outputFolder).Run();


            // TODO: Test output, make sure the values on the results.txt really reflect the number counted in the csv files
            // Parse values from results.txt
            string resultsTextPath = Directory.GetFiles(outputFolder, "allResults.txt", SearchOption.TopDirectoryOnly)
                .FirstOrDefault(); // Try to find the file name "allResults.txt" in the output folder
            if (resultsTextPath is null)
                Assert.Fail("Results file not found.");
            string[] allResultTxtLines = File.ReadAllLines(resultsTextPath); //read all lines from the file
            Assert.That(allResultTxtLines.Length > 0); // make sure there are lines in the file

            //For PSMs
            var allPsmLine = allResultTxtLines.First(p => p.Contains("target PSMs within"));
            int psmCount = int.Parse(allPsmLine.Split(':').Last().Trim());

            //For ProteinGroups
            var proteinGroupLine = allResultTxtLines.First(p => p.Contains("protein groups within"));
            int proteinGroupCount = int.Parse(proteinGroupLine.Split(':').Last().Trim());

            //For GlycoPSMs
            var glycoPsmLine = allResultTxtLines.First(p => p.Contains("O-Glyco PSMs within"));
            int glycoPsmCount = int.Parse(glycoPsmLine.Split(':').Last().Trim()); // read the number of glyco PSMs from the results file

            //For Level1GlycoPSMs
            var level1PsmLine = allResultTxtLines.First(p => p.Contains("Level 1 O-Glyco PSMs within"));
            int level1Psmcount = int.Parse(level1PsmLine.Split(':').Last().Trim()); // read the number of Level1-PSMs from the results file

            // Parse counted number from csv files

            //For PSMs
            var allPsmPath = Path.Combine(outputFolder, "Task", "AllPSMs.psmtsv");
            List<PsmFromTsv> onePercentPsms1 = PsmTsvReader.ReadTsv(allPsmPath, out var errors2)
            .Where(p => p.QValue <= 0.01 && p.DecoyContamTarget != "C" && p.DecoyContamTarget != "D").ToList();
            Assert.That(errors2.Count == 0);// if we cannot find the file, we will get an error message
            int readInPsmsCount = onePercentPsms1.Count;

            //For ProteinGroups
            var allProteinGroupsPath = Path.Combine(outputFolder, "Task", "_AllProteinGroups.tsv");
            string[] proteinGroupHeaders = File.ReadAllLines(allProteinGroupsPath).First().Split("\t");
            int readInProteinCount = File.ReadAllLines(Path.Combine(outputFolder, "Task", "_AllProteinGroups.tsv")).Skip(1)
                .Select(line => line.Split('\t'))
                .Count(p => double.TryParse(p[Array.IndexOf(proteinGroupHeaders, "Protein QValue")], out double qVaule)
                && qVaule < 0.01 && p[Array.IndexOf(proteinGroupHeaders, "Protein Decoy/Contaminant/Target")] != "C"
                && p[Array.IndexOf(proteinGroupHeaders, "Protein Decoy/Contaminant/Target")] != "D");

            //For GlycoPSMs
            string oGlycoPath = Path.Combine(outputFolder, "Task", "oglyco.psmtsv");
            List<PsmFromTsv> onePercentoGlycoPsms = PsmTsvReader.ReadTsv(oGlycoPath, out var errors) //load the PSMs data from the "csv file" and bulid the objects
                .Where(p => p.QValue <= 0.01 && p.DecoyContamTarget != "C" && p.DecoyContamTarget != "D").ToList(); // the filtering (Q<0.01)
            int readInGlycoPsmCount = onePercentoGlycoPsms.Count; // the gPSMs number with Fdr<0.01
            Assert.That(errors.Count == 0);// if we cannot find the file, we will get an error message

            //For Level1GlycoPSMs
            int readInLevel1GlycoPsmCount = onePercentoGlycoPsms.Count(p => p.GlycanLocalizationLevel == EngineLayer.GlycoSearch.LocalizationLevel.Level1); //the level1 gPSMs number

            //Compare the numbers
            Assert.That(psmCount, Is.EqualTo(readInPsmsCount));
            Assert.That(proteinGroupCount, Is.EqualTo(readInProteinCount));
            Assert.That(glycoPsmCount, Is.EqualTo(readInGlycoPsmCount));
            Assert.That(level1Psmcount, Is.EqualTo(readInLevel1GlycoPsmCount));


            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void OGlycoTest_Run5_WriteDecoys() // Test writing decoys, and make sure we can filter the decoys PSMs
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSnip.toml"), MetaMorpheusTask.tomlConfig);
            glycoSearchTask._glycoSearchParameters.WriteContaminants = true;
            glycoSearchTask._glycoSearchParameters.WriteDecoys = true;
            glycoSearchTask._glycoSearchParameters.DecoyType = DecoyType.Reverse;
            glycoSearchTask.CommonParameters = new CommonParameters(dissociationType: DissociationType.HCD, trimMsMsPeaks: false, 
                precursorMassTolerance: new PpmTolerance(6), productMassTolerance: new PpmTolerance(10), qValueThreshold: 1, 
                pepQValueThreshold: 1, scoreCutoff: 1);

            DbForTask targetDbForTask = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoProteinFASTA_7proteins.fasta"), false);
            DbForTask dbContaminant = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P02649.fasta"), true);

            List<string> copiedSpectraFiles = new();
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoPepMix_snip.mzML"); // in order to get enough PSMs to test the filtering, we will copy the spectra file 19 times then get one decpys PSMs in the filterPSMs (Fdr <= 0.01, 1 decoys out of 200 target)
            for(int i = 0; i < 19 ; i++)                    
            {
                var copyPath = Path.Combine(TestContext.CurrentContext.TestDirectory, $@"GlycoTestData\Copy{i}GlycoPepMix_snip.mzML");
                if (!File.Exists(copyPath))
                    File.Copy(spectraFile, copyPath);
                copiedSpectraFiles.Add(copyPath);
            }

            new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) },
                new List<string> { spectraFile, copiedSpectraFiles},
                new List<DbForTask> { targetDbForTask, dbContaminant},
                outputFolder).Run();


            // TODO: Test output, make sure the values on the results.txt really reflect the number counted in the csv files
            // Parse values from results.txt
            string resultsTextPath = Directory.GetFiles(outputFolder, "allResults.txt", SearchOption.TopDirectoryOnly)
                .FirstOrDefault(); // Try to find the file name "allResults.txt" in the output folder
            if (resultsTextPath is null)
                Assert.Fail("Results file not found.");
            string[] allResultTxtLines = File.ReadAllLines(resultsTextPath); //read all lines from the file
            Assert.That(allResultTxtLines.Length > 0); // make sure there are lines in the file

            //For PSMs
            var allPsmLine = allResultTxtLines.First(p => p.Contains("target PSMs within"));
            int psmCount = int.Parse(allPsmLine.Split(':').Last().Trim());

            //For ProteinGroups
            var proteinGroupLine = allResultTxtLines.First(p => p.Contains("protein groups within"));
            int proteinGroupCount = int.Parse(proteinGroupLine.Split(':').Last().Trim());

            //For GlycoPSMs
            var glycoPsmLine = allResultTxtLines.First(p => p.Contains("O-Glyco PSMs within"));
            int glycoPsmCount = int.Parse(glycoPsmLine.Split(':').Last().Trim()); // read the number of glyco PSMs from the results file

            //For Level1GlycoPSMs
            var level1PsmLine = allResultTxtLines.First(p => p.Contains("Level 1 O-Glyco PSMs within"));
            int level1Psmcount = int.Parse(level1PsmLine.Split(':').Last().Trim()); // read the number of Level1-PSMs from the results file

            // Parse counted number from csv files

            //For PSMs
            var allPsmPath = Path.Combine(outputFolder, "Task", "AllPSMs.psmtsv");
            List<PsmFromTsv> onePercentPsms1 = PsmTsvReader.ReadTsv(allPsmPath, out var errors2)
            .Where(p => p.QValue <= 0.01 && p.DecoyContamTarget != "C" && p.DecoyContamTarget != "D").ToList();
            Assert.That(errors2.Count == 0);// if we cannot find the file, we will get an error message
            int readInPsmsCount = onePercentPsms1.Count;

            //For ProteinGroups
            var allProteinGroupsPath = Path.Combine(outputFolder, "Task", "_AllProteinGroups.tsv");
            string[] proteinGroupHeaders = File.ReadAllLines(allProteinGroupsPath).First().Split("\t");
            int readInProteinCount = File.ReadAllLines(Path.Combine(outputFolder, "Task", "_AllProteinGroups.tsv")).Skip(1)
                .Select(line => line.Split('\t'))
                .Count(p => double.TryParse(p[Array.IndexOf(proteinGroupHeaders, "Protein QValue")], out double qVaule)
                && qVaule < 0.01 && p[Array.IndexOf(proteinGroupHeaders, "Protein Decoy/Contaminant/Target")] != "C" //filter the contaminants
                && p[Array.IndexOf(proteinGroupHeaders, "Protein Decoy/Contaminant/Target")] != "D"); // filter the decoys

            //For GlycoPSMs
            string oGlycoPath = Path.Combine(outputFolder, "Task", "oglyco.psmtsv");
            List<PsmFromTsv> onePercentoGlycoPsms = PsmTsvReader.ReadTsv(oGlycoPath, out var errors) //load the PSMs data from the "csv file" and bulid the objects
                .Where(p => p.QValue <= 0.01 && p.DecoyContamTarget != "C" && p.DecoyContamTarget != "D").ToList(); // the filtering (Q<0.01, decoy and contaminat)
            int readInGlycoPsmCount = onePercentoGlycoPsms.Count; // the gPSMs number with Fdr<0.01
            Assert.That(errors.Count == 0);// if we cannot find the file, we will get an error message

            //For Level1GlycoPSMs
            int readInLevel1GlycoPsmCount = onePercentoGlycoPsms.Count(p => p.GlycanLocalizationLevel == EngineLayer.GlycoSearch.LocalizationLevel.Level1); //the level1 gPSMs number

            //Compare the numbers
            Assert.That(psmCount, Is.EqualTo(readInPsmsCount));
            Assert.That(proteinGroupCount, Is.EqualTo(readInProteinCount));
            Assert.That(glycoPsmCount, Is.EqualTo(readInGlycoPsmCount));
            Assert.That(level1Psmcount, Is.EqualTo(readInLevel1GlycoPsmCount));


            copiedSpectraFiles.ForEach(p => File.Delete(p));
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void OGlycoTest_Run6()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSnip.toml"), MetaMorpheusTask.tomlConfig);
            glycoSearchTask._glycoSearchParameters.DoParsimony = false;
            glycoSearchTask._glycoSearchParameters.DoQuantification = true;

            DbForTask targetDbForTask = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoProteinFASTA_7proteins.fasta"), false);
            DbForTask contaminDbForTask = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P13987_contaminant.fasta"), false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoPepMix_snip.mzML");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) }, new List<string> { spectraFile }, new List<DbForTask> { targetDbForTask, contaminDbForTask }, outputFolder).Run();

            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void GlycoTestWithContaminants()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSearchTaskconfigOGlycoTest_Run.toml"), MetaMorpheusTask.tomlConfig);

            DbForTask dbTarget = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P16150.fasta"), false);
            string spectraFileTarget = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf");

            DbForTask dbContaminant = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P02649.fasta"), true);
            string spectraFileContaminant = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\181217_Fusion_(LC2)_NewObj_Serum_deSA_Jacalin_HRM_4h_ETD_HCD_DDA_mz(400_1200)_21707.mgf");

            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) }, new List<string> { spectraFileTarget, spectraFileContaminant }, new List<DbForTask> { dbTarget, dbContaminant }, outputFolder).Run();

            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void GlycoTestWithBadExperimentalDesignFile()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\Task1-GlycoSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            glycoSearchTask._glycoSearchParameters.Normalize = true;

            DbForTask db = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\171025_06_protein.fasta"), false);
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\171025_06subset_1.mzML"), Path.Combine(outputFolder, "171025_06subset_1.mzML"), true);
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\171025_06subset_2.mzML"), Path.Combine(outputFolder, "171025_06subset_2.mzML"), true);

            string spectraFile1 = Path.Combine(outputFolder, "171025_06subset_1.mzML");
            string spectraFile2 = Path.Combine(outputFolder, "171025_06subset_2.mzML");

            List<FlashLFQ.SpectraFileInfo> spectraFiles = new List<FlashLFQ.SpectraFileInfo>();

            //These conditions are such that the experimental design file is bad.
            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(spectraFile1, "condition1", 0, 9, 0));
            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(spectraFile2, "condition2", 5, 0, 4));

            ExperimentalDesign.WriteExperimentalDesignToFile(spectraFiles);

            //this is the crux of the unit test. no experimental design file is passed in and the search task with quant should still run.
            string experimentalDesignFilePath = Path.Combine(outputFolder, "ExperimentalDesign.tsv");

            var readIn = ExperimentalDesign.ReadExperimentalDesign(
                experimentalDesignFilePath,
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(),
                out var errors);

            Assert.That(errors.Count == 1);
            Assert.That(errors[0].Contains("Condition \"condition1\" biorep 1 fraction 1 techrep 1 is missing!"));
            Assert.That(readIn.Count == 2);


            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) }, new List<string> { spectraFile1, spectraFile2 }, new List<DbForTask> { db }, outputFolder).Run();

            List<string> expectedOutput = new()
            {
                "ExperimentalDesign.tsv",
                "_AllProteinGroups.tsv",
                "AllPSMs.psmtsv",
                "AutoGeneratedManuscriptProse.txt",
                "oglyco.psmtsv",
                "protein_oglyco_localization.tsv",
                "results.txt",
                "seen_oglyco_localization.tsv"
            };


            List<string> expectedIndividualFileOutput = new()
            {
                "171025_06subset_1_AllProteinGroups.tsv",
                "171025_06subset_1_AllPSMs.psmtsv",
                "171025_06subset_1oglyco.psmtsv",
                "171025_06subset_1protein_oglyco_localization.tsv",
                "171025_06subset_1seen_oglyco_localization.tsv",
            };

            string outputFolderWithTask = Path.Combine(outputFolder, "Task");
            List<string> output = Directory.GetFiles(outputFolderWithTask).Select(f => Path.GetFileName(f)).ToList();
            List<string> outputFolders = Directory.GetDirectories(outputFolderWithTask).ToList();
            List<string> individualOutputFolders = Directory.GetDirectories(outputFolders.FirstOrDefault()).ToList();
            List<string> individualOutputs = Directory.GetFiles(individualOutputFolders.FirstOrDefault()).Select(f => Path.GetFileName(f)).ToList();

            CollectionAssert.AreEquivalent(expectedOutput, output);
            CollectionAssert.AreEquivalent(expectedIndividualFileOutput, individualOutputs);

            string[] allProteinGroups = File.ReadAllLines(Path.Combine(outputFolderWithTask, "_AllProteinGroups.tsv"));
            string[] proteinGroupFields = allProteinGroups[1].Split('\t');

            Assert.AreEqual("Q9GZM5", proteinGroupFields[0]);

            Directory.Delete(outputFolder, true);
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
        public static void OGlycoTest_GetAllPaths()
        {
            int[] modPos = new int[3] { 2, 4, 6 };
            var glycanBox = OGlycanBoxes[19];
            var boxes = GlycanBox.BuildChildOGlycanBoxes(3, glycanBox.ModIds).ToArray();
            LocalizationGraph localizationGraph = new LocalizationGraph(modPos, glycanBox, boxes, -1);

            for (int i = 0; i < modPos.Length; i++)
            {
                for (int j = 0; j < boxes.Length; j++)
                {
                    localizationGraph.array[i][j] = new AdjNode(i, j, modPos[i], boxes[j]);
                    localizationGraph.array[i][j].CummulativeSources = new List<int> { j }; 
                    localizationGraph.array[i][j].maxCost = 1;
                }
            }
            localizationGraph.array[2][5].CummulativeSources = new List<int> {  4, 5 };

            localizationGraph.array[1][4].CummulativeSources = new List<int> { 2, 4 };

            var allPaths = LocalizationGraph.GetAllHighestScorePaths(localizationGraph.array, boxes);

            Assert.That(allPaths.Count == 3);
            Assert.That(Enumerable.SequenceEqual(allPaths[0], new int[3] { 2, 4, 5 }));
            Assert.That(Enumerable.SequenceEqual(allPaths[1], new int[3] { 4, 4, 5 }));
            Assert.That(Enumerable.SequenceEqual(allPaths[2], new int[3] { 5, 5, 5 }));

            var firstPath = LocalizationGraph.GetFirstPath(localizationGraph.array, boxes);  
            Assert.That(Enumerable.SequenceEqual(firstPath, new int[3] { 2, 4, 5 }));
        }

        [Test]
        public static void OGlycoTest_GetLocalizedPath()
        {
            int[] modPos = new int[1] { 4 };
            var glycanBox = OGlycanBoxes[1];
            var boxes = GlycanBox.BuildChildOGlycanBoxes(1, glycanBox.ModIds).ToArray();
            LocalizationGraph localizationGraph = new LocalizationGraph(modPos, glycanBox, boxes, -1);

            for (int i = 0; i < modPos.Length; i++)
            {
                for (int j = 0; j < boxes.Length; j++)
                {
                    localizationGraph.array[i][j] = new AdjNode(i, j, modPos[i], boxes[j]);
                    localizationGraph.array[i][j].CummulativeSources = new List<int> { j };
                }
            }
            localizationGraph.TotalScore = 1;

           var allPaths = LocalizationGraph.GetAllHighestScorePaths(localizationGraph.array, boxes);

            var route = LocalizationGraph.GetLocalizedPath(localizationGraph, allPaths.First());

            Assert.That(route.Mods.First().Item3);
        }

        [Test]
        public static void OGlycoTest_DissociationTypeContainETD()
        {
            DissociationType dissociationType = DissociationType.Custom;

            List<ProductType> customIons = new List<ProductType> { ProductType.zDot };

            var isETDType = GlycoPeptides.DissociationTypeContainETD(dissociationType, customIons);

            Assert.That(isETDType);
        }

        [Test]
        [TestCase(true, true, false, true)] //write decoys and contaminants, database is targets, adds protein stuff
        [TestCase(true, true, true, false)] //write decoys and contaminants, database is contaminants, all contaminants so protein stuff is skipped
        [TestCase(false, false, false, true)] //do NOT write decoys and contaminants, database is targets, adds protein stuff
        public static void OGlycoTest_ProteinInference(bool writeContaminants, bool writeDecoys, bool isContaminant, bool addProteinOutput)
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            string proteinDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\N_O_glycoWithFileSpecific\\FourMucins_NoSigPeps_FASTA.fasta");
            string spectraFileDirctory = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\N_O_glycoWithFileSpecific");
            List<string> rawFilePaths = Directory.GetFiles(spectraFileDirctory).Where(p => p.Contains("mzML")).ToList();

            // run task
            CommonParameters commonParameters = new(dissociationType: DissociationType.HCD, ms2childScanDissociationType: DissociationType.EThcD);

            Directory.CreateDirectory(outputFolder);
            var glycoSearchTask = new GlycoSearchTask()
            {
                CommonParameters = commonParameters,
                _glycoSearchParameters = new GlycoSearchParameters()
                {
                    OGlycanDatabasefile = "OGlycan.gdb",
                    NGlycanDatabasefile = "NGlycan.gdb",
                    GlycoSearchType = GlycoSearchType.OGlycanSearch,
                    OxoniumIonFilt = true,
                    DecoyType = DecoyType.Reverse,
                    GlycoSearchTopNum = 50,
                    MaximumOGlycanAllowed = 4,
                    DoParsimony = true,
                    WriteContaminants = writeContaminants,
                    WriteDecoys = writeDecoys
                }
            };
            glycoSearchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(proteinDatabase, isContaminant) }, rawFilePaths, "");

            List<string> expectedOutput = new()
            {
                "_AllProteinGroups.tsv",
                "AutoGeneratedManuscriptProse.txt",
                "results.txt",
                "oglyco.psmtsv",
                "AllPSMs.psmtsv"
            };
            if (addProteinOutput)
            {
                expectedOutput.Add("protein_oglyco_localization.tsv");
                expectedOutput.Add("seen_oglyco_localization.tsv");
            }

            List<string> output = Directory.GetFiles(outputFolder).Select(f => Path.GetFileName(f)).ToList();

            CollectionAssert.IsSubsetOf(expectedOutput, output);

            string[] allProteinGroups = File.ReadAllLines(Path.Combine(outputFolder, "_AllProteinGroups.tsv"));
            string[] proteinGroupFields = allProteinGroups[1].Split('\t');

            Assert.AreEqual("Q8WXI7.3", proteinGroupFields[0]);

            Directory.Delete(outputFolder, true);
        }
        [Test]
        public static void GlycoSearchIndividualFileFolderOutputTest()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            string proteinDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\N_O_glycoWithFileSpecific\\FourMucins_NoSigPeps_FASTA.fasta");
            string spectraFileDirctory = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\N_O_glycoWithFileSpecific");
            List<string> rawFilePaths = Directory.GetFiles(spectraFileDirctory).Where(p => p.Contains("mzML")).ToList();

            // run task
            CommonParameters commonParameters = new(dissociationType: DissociationType.HCD, ms2childScanDissociationType: DissociationType.EThcD);

            Directory.CreateDirectory(outputFolder);
            var glycoSearchTask = new GlycoSearchTask()
            {
                CommonParameters = commonParameters,
                _glycoSearchParameters = new GlycoSearchParameters()
                {
                    OGlycanDatabasefile = "OGlycan.gdb",
                    NGlycanDatabasefile = "NGlycan.gdb",
                    GlycoSearchType = GlycoSearchType.OGlycanSearch,
                    OxoniumIonFilt = true,
                    DecoyType = DecoyType.Reverse,
                    GlycoSearchTopNum = 50,
                    MaximumOGlycanAllowed = 4,
                    DoParsimony = true,
                    WriteContaminants = true,
                    WriteDecoys = true,
                    WriteIndividualFiles = true
                }
            };
            glycoSearchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(proteinDatabase, false) }, rawFilePaths, "");

            List<string> expectedOutput = new()
            {
                "AutoGeneratedManuscriptProse.txt",
                "results.txt",
                "oglyco.psmtsv",
                "AllPSMs.psmtsv",
                "_AllProteinGroups.tsv",
                "protein_oglyco_localization.tsv",
                "seen_oglyco_localization.tsv"
            };

            List<string> expectedIndividualFileOutput = new()
            {
                "2019_09_16_StcEmix_35trig_EThcD25_rep1_4999-5968_AllProteinGroups.tsv",
                "2019_09_16_StcEmix_35trig_EThcD25_rep1_4999-5968_AllPSMs.psmtsv",
                "2019_09_16_StcEmix_35trig_EThcD25_rep1_4999-5968oglyco.psmtsv",
                "2019_09_16_StcEmix_35trig_EThcD25_rep1_4999-5968protein_oglyco_localization.tsv",
                "2019_09_16_StcEmix_35trig_EThcD25_rep1_4999-5968seen_oglyco_localization.tsv"
            };



            List<string> output = Directory.GetFiles(outputFolder).Select(f => Path.GetFileName(f)).ToList();
            List<string> outputFolders = Directory.GetDirectories(outputFolder).ToList();
            List<string> individualOutputFolders = Directory.GetDirectories(outputFolders.FirstOrDefault()).ToList();
            List<string> individualOutputs = Directory.GetFiles(individualOutputFolders.FirstOrDefault()).Select(f => Path.GetFileName(f)).ToList();

            CollectionAssert.IsSubsetOf(expectedOutput, output);
            CollectionAssert.IsSubsetOf(expectedIndividualFileOutput, individualOutputs);

            string[] allProteinGroups = File.ReadAllLines(Path.Combine(outputFolder, "_AllProteinGroups.tsv"));
            string[] proteinGroupFields = allProteinGroups[1].Split('\t');

            Assert.AreEqual("Q8WXI7.3", proteinGroupFields[0]);

            Directory.Delete(outputFolder, true);
        }
        [Test]
        public static void NandO_GlycoSearchIndividualFileFolderOutputTest()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            string proteinDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\N_O_glycoWithFileSpecific\\FourMucins_NoSigPeps_FASTA.fasta");
            string spectraFileDirctory = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\N_O_glycoWithFileSpecific");
            List<string> rawFilePaths = Directory.GetFiles(spectraFileDirctory).Where(p => p.Contains("mzML")).ToList();

            // run task
            CommonParameters commonParameters = new(dissociationType: DissociationType.HCD, ms2childScanDissociationType: DissociationType.EThcD);

            Directory.CreateDirectory(outputFolder);
            var glycoSearchTask = new GlycoSearchTask()
            {
                CommonParameters = commonParameters,
                _glycoSearchParameters = new GlycoSearchParameters()
                {
                    OGlycanDatabasefile = "OGlycan.gdb",
                    NGlycanDatabasefile = "NGlycan.gdb",
                    GlycoSearchType = GlycoSearchType.N_O_GlycanSearch,
                    OxoniumIonFilt = true,
                    DecoyType = DecoyType.Reverse,
                    GlycoSearchTopNum = 50,
                    MaximumOGlycanAllowed = 4,
                    DoParsimony = true,
                    WriteContaminants = true,
                    WriteDecoys = true,
                    WriteIndividualFiles = true
                }
            };
            glycoSearchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(proteinDatabase, false) }, rawFilePaths, "");

            List<string> expectedOutput = new()
            {
                "_AllProteinGroups.tsv",
                "AutoGeneratedManuscriptProse.txt",
                "results.txt",
                "no_glyco.psmtsv",
                "AllPSMs.psmtsv",
                "protein_no_glyco_localization.tsv",
                "seen_no_glyco_localization.tsv"
            };

            List<string> expectedIndividualFileOutput = new()
            {
                "2019_09_16_StcEmix_35trig_EThcD25_rep1_4999-5968_AllProteinGroups.tsv",
                "2019_09_16_StcEmix_35trig_EThcD25_rep1_4999-5968_AllPSMs.psmtsv",
                "2019_09_16_StcEmix_35trig_EThcD25_rep1_4999-5968no_glyco.psmtsv",
                "2019_09_16_StcEmix_35trig_EThcD25_rep1_4999-5968protein_no_glyco_localization.tsv",
                "2019_09_16_StcEmix_35trig_EThcD25_rep1_4999-5968seen_no_glyco_localization.tsv"
            };



            List<string> output = Directory.GetFiles(outputFolder).Select(f => Path.GetFileName(f)).ToList();
            List<string> outputFolders = Directory.GetDirectories(outputFolder).ToList();
            List<string> individualOutputFolders = Directory.GetDirectories(outputFolders.FirstOrDefault()).ToList();
            List<string> individualOutputs = Directory.GetFiles(individualOutputFolders.FirstOrDefault()).Select(f => Path.GetFileName(f)).ToList();

            CollectionAssert.IsSubsetOf(expectedOutput, output);
            CollectionAssert.IsSubsetOf(expectedIndividualFileOutput, individualOutputs);

            string[] allProteinGroups = File.ReadAllLines(Path.Combine(outputFolder, "_AllProteinGroups.tsv"));
            string[] proteinGroupFields = allProteinGroups[1].Split('\t');

            Assert.AreEqual("Q8WXI7.3", proteinGroupFields[0]);

            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void N_GlycoSearchIndividualFileFolderOutputTest()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            string proteinDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\N_O_glycoWithFileSpecific\\FourMucins_NoSigPeps_FASTA.fasta");
            string spectraFileDirctory = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\N_O_glycoWithFileSpecific");
            List<string> rawFilePaths = Directory.GetFiles(spectraFileDirctory).Where(p => p.Contains("mzML")).ToList();

            // run task
            CommonParameters commonParameters = new(dissociationType: DissociationType.HCD, ms2childScanDissociationType: DissociationType.EThcD);

            Directory.CreateDirectory(outputFolder);
            var glycoSearchTask = new GlycoSearchTask()
            {
                CommonParameters = commonParameters,
                _glycoSearchParameters = new GlycoSearchParameters()
                {
                    OGlycanDatabasefile = "OGlycan.gdb",
                    NGlycanDatabasefile = "NGlycan.gdb",
                    GlycoSearchType = GlycoSearchType.NGlycanSearch,
                    OxoniumIonFilt = true,
                    DecoyType = DecoyType.Reverse,
                    GlycoSearchTopNum = 50,
                    MaximumOGlycanAllowed = 4,
                    DoParsimony = true,
                    WriteContaminants = true,
                    WriteDecoys = true,
                    WriteIndividualFiles = true,
                    DoQuantification = true
                }
            };
            glycoSearchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(proteinDatabase, false) }, rawFilePaths, "");

            List<string> expectedOutput = new()
            {
                "_AllProteinGroups.tsv",
                "AllQuantifiedPeaks.tsv",
                "AllQuantifiedPeptides.tsv",
                "AllQuantifiedProteins.tsv",
                "AutoGeneratedManuscriptProse.txt",
                "results.txt",
                "ExperimentalDesign.tsv",
                "nglyco.psmtsv",
                "AllPSMs.psmtsv",
            };

            List<string> expectedIndividualFileOutput = new()
            {
                "2019_09_16_StcEmix_35trig_EThcD25_rep1_4999-5968nglyco.psmtsv",
                "2019_09_16_StcEmix_35trig_EThcD25_rep1_4999-5968_AllProteinGroups.tsv",
                "2019_09_16_StcEmix_35trig_EThcD25_rep1_4999-5968_AllPSMs.psmtsv",
                "2019_09_16_StcEmix_35trig_EThcD25_rep1_4999-5968_QuantifiedPeaks.tsv",
                "2019_09_16_StcEmix_35trig_EThcD25_rep1_4999-5968_QuantifiedPeptides.tsv",
                "2019_09_16_StcEmix_35trig_EThcD25_rep1_4999-5968_QuantifiedProteins.tsv"
            };



            List<string> output = Directory.GetFiles(outputFolder).Select(f => Path.GetFileName(f)).ToList();
            List<string> outputFolders = Directory.GetDirectories(outputFolder).ToList();
            List<string> individualOutputFolders = Directory.GetDirectories(outputFolders.FirstOrDefault()).ToList();
            List<string> individualOutputs = Directory.GetFiles(individualOutputFolders.FirstOrDefault()).Select(f => Path.GetFileName(f)).ToList();

            CollectionAssert.IsSubsetOf(expectedOutput, output);
            CollectionAssert.IsSubsetOf(expectedIndividualFileOutput, individualOutputs);

            string[] allProteinGroups = File.ReadAllLines(Path.Combine(outputFolder, "_AllProteinGroups.tsv"));
            string[] proteinGroupFields = allProteinGroups[1].Split('\t');

            Assert.AreEqual("Q8WXI7.3", proteinGroupFields[0]);

            Directory.Delete(outputFolder, true);
        }

        //Test oglycopeptide quant with two conditions, experimental design and normalization
        [Test]
        public static void TestGlycoQuant()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\Task1-GlycoSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            glycoSearchTask._glycoSearchParameters.Normalize = true;

            DbForTask db = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\171025_06_protein.fasta"), false);
            string spectraFile1 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\171025_06subset_1.mzML");
            string spectraFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\171025_06subset_2.mzML");

            List<FlashLFQ.SpectraFileInfo> spectraFiles = new List<FlashLFQ.SpectraFileInfo>();

            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(spectraFile1, "condition1", 0, 0, 0));
            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(spectraFile2, "condition2", 0, 0, 0));

            ExperimentalDesign.WriteExperimentalDesignToFile(spectraFiles);

            string experimentalDesignFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\ExperimentalDesign.tsv");

            var readIn = ExperimentalDesign.ReadExperimentalDesign(
                experimentalDesignFilePath,
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(),
                out var errors);

            Assert.That(!errors.Any());
            Assert.That(readIn.Count == 2);


            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) }, new List<string> { spectraFile1, spectraFile2 }, new List<DbForTask> { db }, outputFolder).Run();

            List<string> expectedOutput = new()
            {
                "ExperimentalDesign.tsv",
                "_AllProteinGroups.tsv",
                "AllPSMs.psmtsv",
                "AllQuantifiedPeaks.tsv",
                "AllQuantifiedPeptides.tsv",
                "AllQuantifiedProteins.tsv",
                "AutoGeneratedManuscriptProse.txt",
                "oglyco.psmtsv",
                "protein_oglyco_localization.tsv",
                "results.txt",
                "seen_oglyco_localization.tsv"
            };


            List<string> expectedIndividualFileOutput = new()
            {
                "171025_06subset_1_AllProteinGroups.tsv",
                "171025_06subset_1_AllPSMs.psmtsv",
                "171025_06subset_1_QuantifiedPeaks.tsv",
                "171025_06subset_1_QuantifiedPeptides.tsv",
                "171025_06subset_1_QuantifiedProteins.tsv",
                "171025_06subset_1oglyco.psmtsv",
                "171025_06subset_1protein_oglyco_localization.tsv",
                "171025_06subset_1seen_oglyco_localization.tsv",
            };


            string outputFolderWithTask = Path.Combine(outputFolder, "Task");
            List<string> output = Directory.GetFiles(outputFolderWithTask).Select(f => Path.GetFileName(f)).ToList();
            List<string> outputFolders = Directory.GetDirectories(outputFolderWithTask).ToList();
            List<string> individualOutputFolders = Directory.GetDirectories(outputFolders.FirstOrDefault()).ToList();
            List<string> individualOutputs = Directory.GetFiles(individualOutputFolders.FirstOrDefault()).Select(f => Path.GetFileName(f)).ToList();

            CollectionAssert.AreEquivalent(expectedOutput, output);
            CollectionAssert.AreEquivalent(expectedIndividualFileOutput, individualOutputs);

            string[] allProteinGroups = File.ReadAllLines(Path.Combine(outputFolderWithTask, "_AllProteinGroups.tsv"));
            string[] proteinGroupFields = allProteinGroups[1].Split('\t');

            Assert.AreEqual("Q9GZM5", proteinGroupFields[0]);

            File.Delete(experimentalDesignFilePath);
            Directory.Delete(outputFolder, true);
        }
        //Test oglycopeptide quant with two conditions, experimental design and NO-normalization and NO protein inference
        [Test]
        public static void TestGlycoQuant2()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\Task1-GlycoSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            glycoSearchTask._glycoSearchParameters.DoParsimony = false;

            DbForTask db = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\171025_06_protein.fasta"), false);
            string spectraFile1 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\171025_06subset_1.mzML");
            string spectraFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\171025_06subset_2.mzML");

            List<FlashLFQ.SpectraFileInfo> spectraFiles = new List<FlashLFQ.SpectraFileInfo>();

            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(spectraFile1, "condition1", 0, 0, 0));
            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(spectraFile2, "condition2", 0, 0, 0));

            ExperimentalDesign.WriteExperimentalDesignToFile(spectraFiles);

            string experimentalDesignFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\ExperimentalDesign.tsv");

            var readIn = ExperimentalDesign.ReadExperimentalDesign(
                experimentalDesignFilePath,
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(),
                out var errors);

            Assert.That(!errors.Any());
            Assert.That(readIn.Count == 2);


            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) }, new List<string> { spectraFile1, spectraFile2 }, new List<DbForTask> { db }, outputFolder).Run();

            List<string> expectedOutput = new()
            {
                "ExperimentalDesign.tsv",
                "AllPSMs.psmtsv",
                "AllQuantifiedPeaks.tsv",
                "AllQuantifiedPeptides.tsv",
                "AllQuantifiedProteins.tsv",
                "AutoGeneratedManuscriptProse.txt",
                "oglyco.psmtsv",
                "protein_oglyco_localization.tsv",
                "results.txt",
                "seen_oglyco_localization.tsv"
            };


            List<string> expectedIndividualFileOutput = new()
            {
                "171025_06subset_1_AllPSMs.psmtsv",
                "171025_06subset_1_QuantifiedPeaks.tsv",
                "171025_06subset_1_QuantifiedPeptides.tsv",
                "171025_06subset_1_QuantifiedProteins.tsv",
                "171025_06subset_1oglyco.psmtsv",
                "171025_06subset_1protein_oglyco_localization.tsv",
                "171025_06subset_1seen_oglyco_localization.tsv",
            };


            string outputFolderWithTask = Path.Combine(outputFolder, "Task");
            List<string> output = Directory.GetFiles(outputFolderWithTask).Select(f => Path.GetFileName(f)).ToList();
            List<string> outputFolders = Directory.GetDirectories(outputFolderWithTask).ToList();
            List<string> individualOutputFolders = Directory.GetDirectories(outputFolders.FirstOrDefault()).ToList();
            List<string> individualOutputs = Directory.GetFiles(individualOutputFolders.FirstOrDefault()).Select(f => Path.GetFileName(f)).ToList();

            CollectionAssert.AreEquivalent(expectedOutput, output);
            CollectionAssert.AreEquivalent(expectedIndividualFileOutput, individualOutputs);

            File.Delete(experimentalDesignFilePath);
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void GlycoQuantWithNoExperimentalDesignFileTest()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\Task1-GlycoSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            glycoSearchTask._glycoSearchParameters.Normalize = true;

            DbForTask db = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\171025_06_protein.fasta"), false);
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\171025_06subset_1.mzML"), Path.Combine(outputFolder,"171025_06subset_1.mzML"), true);
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\171025_06subset_2.mzML"), Path.Combine(outputFolder, "171025_06subset_2.mzML"), true);


            string spectraFile1 = Path.Combine(outputFolder, "171025_06subset_1.mzML");
            string spectraFile2 = Path.Combine(outputFolder, "171025_06subset_2.mzML");

            List<FlashLFQ.SpectraFileInfo> spectraFiles = new List<FlashLFQ.SpectraFileInfo>();

            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(spectraFile1, "condition1", 0, 0, 0));
            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(spectraFile2, "condition2", 0, 0, 0));

            ExperimentalDesign.WriteExperimentalDesignToFile(spectraFiles);

            //this is the crux of the unit test. no experimental design file is passed in and the search task with quant should still run.
            string experimentalDesignFilePath = null;

            var readIn = ExperimentalDesign.ReadExperimentalDesign(
                experimentalDesignFilePath,
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(),
                out var errors);

            Assert.That(errors.Count == 1);
            Assert.That(errors[0].Contains("Experimental design file not found"));
            Assert.That(readIn.Count == 0);


            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) }, new List<string> { spectraFile1, spectraFile2 }, new List<DbForTask> { db }, outputFolder).Run();

            List<string> expectedOutput = new()
            {
                "ExperimentalDesign.tsv",
                "AllQuantifiedProteins.tsv",
                "AllPSMs.psmtsv",
                "AllQuantifiedPeaks.tsv",
                "AllQuantifiedPeptides.tsv",
                "AutoGeneratedManuscriptProse.txt",
                "oglyco.psmtsv",
                "protein_oglyco_localization.tsv",
                "results.txt",
                "seen_oglyco_localization.tsv",
                "_AllProteinGroups.tsv"
            };


            List<string> expectedIndividualFileOutput = new()
            {
                "171025_06subset_1_AllProteinGroups.tsv",
                "171025_06subset_1_AllPSMs.psmtsv",
                "171025_06subset_1_QuantifiedPeaks.tsv",
                "171025_06subset_1_QuantifiedPeptides.tsv",
                "171025_06subset_1_QuantifiedProteins.tsv",
                "171025_06subset_1oglyco.psmtsv",
                "171025_06subset_1protein_oglyco_localization.tsv",
                "171025_06subset_1seen_oglyco_localization.tsv",
            };

            string outputFolderWithTask = Path.Combine(outputFolder, "Task");
            List<string> output = Directory.GetFiles(outputFolderWithTask).Select(f => Path.GetFileName(f)).ToList();
            List<string> outputFolders = Directory.GetDirectories(outputFolderWithTask).ToList();
            List<string> individualOutputFolders = Directory.GetDirectories(outputFolders.FirstOrDefault()).ToList();
            List<string> individualOutputs = Directory.GetFiles(individualOutputFolders.FirstOrDefault()).Select(f => Path.GetFileName(f)).ToList();

            CollectionAssert.AreEquivalent(expectedOutput, output);
            CollectionAssert.AreEquivalent(expectedIndividualFileOutput, individualOutputs);

            string[] allProteinGroups = File.ReadAllLines(Path.Combine(outputFolderWithTask, "AllQuantifiedProteins.tsv"));
            string[] proteinGroupFields = allProteinGroups[2].Split('\t');

            Assert.AreEqual("Q9GZM5", proteinGroupFields[0]);

            Directory.Delete(outputFolder, true);
        }

        //Test oglycopeptide quant with two conditions, experimental design and normalization
        [Test]
        public static void TestExperimentalDesignError()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\Task1-GlycoSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            glycoSearchTask._glycoSearchParameters.Normalize = true;

            DbForTask db = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\171025_06_protein.fasta"), false);
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\171025_06subset_1.mzML"), Path.Combine(outputFolder, "171025_06subset_1.mzML"), true);
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\171025_06subset_2.mzML"), Path.Combine(outputFolder, "171025_06subset_2.mzML"), true);


            string spectraFile1 = Path.Combine(outputFolder, "171025_06subset_1.mzML");
            string spectraFile2 = Path.Combine(outputFolder, "171025_06subset_2.mzML");

            CommonParameters commonParameters = new CommonParameters(
                taskDescriptor: "AveragingTask",
                dissociationType: DissociationType.LowCID,
                maxThreadsToUsePerFile: 2);

            SpectralAveragingParameters options = new SpectralAveragingParameters()
            {
                SpectraFileAveragingType = SpectraFileAveragingType.AverageDdaScans,
                NumberOfScansToAverage = 5,
            };
            SpectralAveragingTask averagingTask = new(options) { CommonParameters = commonParameters };

            // run the tasks
            EverythingRunnerEngine everythingRunnerEngine = new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("Task1-Average", averagingTask), ("Task2-Search", glycoSearchTask) },
                new List<string> { spectraFile1, spectraFile2 },
                new List<DbForTask> { db },
                outputFolder);

            // set up original experimental design (input to calibration)
            SpectraFileInfo fileInfo1 = new SpectraFileInfo(spectraFile1, "condition", 0, 0, 0);
            SpectraFileInfo fileInfo2 = new SpectraFileInfo(spectraFile2, "condition", 0, 0, 0);

            var experimentalDesignFilePath = ExperimentalDesign.WriteExperimentalDesignToFile(new List<SpectraFileInfo> { fileInfo1, fileInfo2 });

            using (StreamWriter writer = new StreamWriter(File.OpenWrite(experimentalDesignFilePath)))
            {
                writer.WriteLine("\t2\t4\t7");
            }

            everythingRunnerEngine.Run();

            // make sure everything is there
            Assert.That(File.Exists(Path.Combine(outputFolder, "allResults.txt")));

            Assert.That(Directory.Exists(Path.Combine(outputFolder, "Task Settings")));
            Assert.That(File.Exists(Path.Combine(outputFolder, "Task Settings", "Task1-Averageconfig.toml")));
            Assert.That(File.Exists(Path.Combine(outputFolder, "Task Settings", "Task2-Searchconfig.toml")));

            Assert.That(Directory.Exists(Path.Combine(outputFolder, "Task1-Average")));
            Assert.That(File.Exists(Path.Combine(outputFolder, "Task1-Average", "AutoGeneratedManuscriptProse.txt")));
            Assert.That(File.Exists(Path.Combine(outputFolder, "Task1-Average", "results.txt")));
            Assert.That(File.Exists(Path.Combine(outputFolder, "Task1-Average", "171025_06subset_1-averaged.mzML")));
            Assert.That(File.Exists(Path.Combine(outputFolder, "Task1-Average", "171025_06subset_2-averaged.mzML")));

            Assert.That(Directory.Exists(Path.Combine(outputFolder, "Task2-Search")));
            Assert.That(File.Exists(Path.Combine(outputFolder, "Task2-Search", "AutoGeneratedManuscriptProse.txt")));
            Assert.That(File.Exists(Path.Combine(outputFolder, "Task2-Search", "results.txt")));
            Assert.That(File.Exists(Path.Combine(outputFolder, "Task2-Search", "AllPSMs.psmtsv")));
            Assert.That(File.Exists(Path.Combine(outputFolder, "Task2-Search", "AllQuantifiedPeaks.tsv")));
            Assert.That(File.Exists(Path.Combine(outputFolder, "Task2-Search", "AllQuantifiedPeptides.tsv")));
            Assert.That(File.Exists(Path.Combine(outputFolder, "Task2-Search", "oglyco.psmtsv")));
            Assert.That(File.Exists(Path.Combine(outputFolder, "Task2-Search", "seen_oglyco_localization.tsv")));
            Assert.That(File.Exists(Path.Combine(outputFolder, "Task2-Search", "protein_oglyco_localization.tsv")));

            File.Delete(experimentalDesignFilePath);
            Directory.Delete(outputFolder, true);
        }
        [Test]
        [TestCase(false, 2, 1, 1)] // pre output: 1 intensity column, post output: 2 intensity column
        [TestCase(true, 2, 3, 1)]  // pre output: 1 intensity column, post output: 6 intensity column
        [TestCase(true, 2, 3, 2)]  // pre output: 1 intensity column, post output: 12 intensity column


        public static void TestGlycoProteinQuantFileHeaders(bool hasDefinedExperimentalDesign, int bioreps, int fractions, int techreps)
        {
            string condition = hasDefinedExperimentalDesign ? "TestCondition" : "";
            List<SpectraFileInfo> fileInfos = new();

            // create the unit test directory
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            // create the .mzML files to search/quantify
            for (int b = 0; b < bioreps; b++)
            {
                for (int f = 0; f < fractions; f++)
                {
                    for (int r = 0; r < techreps; r++)
                    {
                        string fileToWrite = "file_" + "b" + b + "f" + f + "r" + r + ".mzML";
                        string fullPath = Path.Combine(outputFolder, fileToWrite);
                        SpectraFileInfo spectraFileInfo = new(fullPath, condition, b, r, f);

                        File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\171025_06subset_1.mzML"), fullPath, true);


                        fileInfos.Add(spectraFileInfo);
                    }
                }
            }

            // write the experimental design for this quantification test
            if (hasDefinedExperimentalDesign)
            {
                _ = ExperimentalDesign.WriteExperimentalDesignToFile(fileInfos);
            }

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\Task1-GlycoSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            glycoSearchTask._glycoSearchParameters.DoParsimony = true;

            DbForTask db = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\QuantData\171025_06_protein.fasta"), false);

            glycoSearchTask.RunTask(outputFolder, new List<DbForTask> { db }, fileInfos.Select(p => p.FullFilePathWithExtension).ToList(), "");

            // read in the protein quant results
            Assert.That(File.Exists(Path.Combine(outputFolder, "AllQuantifiedPeptides.tsv")));
            string[] lines = File.ReadAllLines(Path.Combine(outputFolder, "AllQuantifiedPeptides.tsv"));

            // check the intensity column headers
            List<string> splitHeader = lines[0].Split(new char[] { '\t' }).ToList();
            List<string> intensityColumnHeaders = splitHeader.Where(p => p.Contains("Intensity", StringComparison.OrdinalIgnoreCase)).ToList();

            Assert.That(intensityColumnHeaders.Count == bioreps* fractions* techreps); // We change the search funtion allowed to get the PSMs from the duplicate file.
                                                                                       // Ex. we have 2 bioreps, 3 fractions, 1 techrep, then we get 6 intensity columns
            Directory.Delete(outputFolder, true);
        }
        [Test]
        public static void TestRunSpecificPostGlycoSearchAnalysis()
        {
            //code coverage unit test for an unused abstract method in post search analysis
            var task = new PostGlycoSearchAnalysisTask();
            task.RunTask(TestContext.CurrentContext.TestDirectory, new List<DbForTask>(), new List<string>(), "");
        }

        [Test]
        public static void TestWriteSpectrumLibrary()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSearchTaskconfigOGlycoTest_Run.toml"), MetaMorpheusTask.tomlConfig);
            glycoSearchTask._glycoSearchParameters.WriteSpectrumLibrary =  true;

            DbForTask db = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P16150.fasta"), false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) }, new List<string> { spectraFile }, new List<DbForTask> { db }, outputFolder).Run();

            List<string> folderContents = Directory.GetFiles(Path.Combine(outputFolder,"Task")).Select(p => Path.GetFileName(p)).ToList();
            string spectrumLibrary = folderContents.Where(p => p.Contains("SpectralLibrary")).FirstOrDefault();
            Assert.That(spectrumLibrary.Contains(".msp"));

            Directory.Delete(outputFolder, true);
        }
    }
}
