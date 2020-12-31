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
using EngineLayer.GlycoSearch;
using MathNet.Numerics.LinearRegression;

namespace Test
{
    [TestFixture]
    public class XLTestOGlyco
    {
        private static GlycanBox[] OGlycanBoxes { get; set; }

        [OneTimeSetUp]
        public static void Setup()
        {
            GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.OGlycanLocations.Where(p => p.Contains("OGlycan.gdb")).First(), true, true).ToArray();
            GlycanBox.GlobalOGlycanMods = GlycanBox.BuildGlobalOGlycanMods(GlycanBox.GlobalOGlycans).ToArray();
            OGlycanBoxes = GlycanBox.BuildGlycanBoxes(3, GlycanBox.GlobalOGlycans, GlycanBox.GlobalOGlycanMods).OrderBy(p => p.Mass).ToArray();
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
        public static void GlycoTest_MotifExist()
        {
            string baseSeq = "AFGQFFSPGEVIYNKTDRAG";
            var exist = GlycoPeptides.MotifExist(baseSeq, new string[] { "Nxt", "Nxs" });
            Assert.That(exist);
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
            //This test is to test ETD on proline. As proline didn't have proline related product ions.

            //Get glycanBox
            var glycanBox = OGlycanBoxes[8];

            Protein protein = new Protein("PTLFKNVSLYK", "");
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();

            List<int> modPos = GlycoPeptides.GetPossibleModSites(peptide, new string[] { "S", "T" });

            var peptideWithMod = GlycoPeptides.GlyGetTheoreticalPeptide(modPos.ToArray(), peptide, OGlycanBoxes[8], GlycanBox.GlobalOGlycanMods);
            Assert.That(peptideWithMod.FullSequence == "PT[O-Glycosylation:N1A1 on X]LFKNVS[O-Glycosylation:N1 on X]LYK");

            var fragments_hcd = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.HCD, peptide, peptideWithMod);

            var fragments_ethcd = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.EThcD, peptide, peptideWithMod);
        }

        [Test]
        public static void OGlycoTest_FragmentIons2()
        {
            //Get glycanBox
            var glycanBox = OGlycanBoxes[24];

            Protein protein = new Protein("TVYLGASK", "");
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();

            List<int> modPos = new List<int> { 2, 8 };

            var peptideWithMod = GlycoPeptides.GlyGetTheoreticalPeptide(modPos.ToArray(), peptide, OGlycanBoxes[24], GlycanBox.GlobalOGlycanMods);
            Assert.That(peptideWithMod.FullSequence == "T[O-Glycosylation:H1N1 on X]VYLGAS[O-Glycosylation:H1N1A1 on X]K");

            var fragments_etd = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.ETD, peptide, peptideWithMod);


            Assert.That(fragments_etd.Count == 22);
            Assert.That(fragments_etd.Last().Annotation == "zDot8");
            Assert.That(fragments_etd.Last().NeutralMass > 1824);

            var fragments_hcd = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.HCD, peptide, peptideWithMod);
            Assert.That(fragments_hcd.Where(p=>p.ProductType == ProductType.M).Count() == 5);

            var fragments_ethcd = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.EThcD, peptide, peptideWithMod);
            Assert.That(fragments_ethcd.Where(p => p.ProductType == ProductType.M).Count() == 5);
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

            int[] modPos = GlycoPeptides.GetPossibleModSites(peptide, new string[] { "S", "T" }).OrderBy(v=>v).ToArray();
            var boxes = GlycanBox.BuildChildGlycanBoxes(glycanBox.ModIds, GlycanBox.GlobalOGlycans, GlycanBox.GlobalOGlycanMods).ToArray();
            Assert.That(boxes.Count() == 6);

            //Get Unlocal Fragment
            var unlocalCost = GlycoPeptides.GetUnlocalFragment(products, modPos, glycanBox);
            Assert.That(unlocalCost.Count == 5); //Basicly, the unlocal are c/z ions that don't localize glycosylation. 

            //Get scan
            CommonParameters commonParameters = new CommonParameters(dissociationType: DissociationType.EThcD, trimMsMsPeaks:false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_4565.mgf");
            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
            var scans = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).ToArray();

            //Known peptideWithMod match.
            var peptideWithMod = GlycoPeptides.GlyGetTheoreticalPeptide(new int[3] { 10, 2, 3}, peptide, glycanBox, GlycanBox.GlobalOGlycanMods);
            Assert.That(peptideWithMod.FullSequence == "T[O-Glycosylation:H1N1 on X]T[O-Glycosylation:H1N1 on X]GSLEPSS[O-Glycosylation:N1 on X]GASGPQVSSVK");
            List<Product> knownProducts = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.EThcD, peptide, peptideWithMod);
            var matchedKnownFragmentIons = MetaMorpheusEngine.MatchFragmentIons(scans.First(), knownProducts, commonParameters);

            //Graph Localization
            LocalizationGraph localizationGraph = new LocalizationGraph(modPos, glycanBox, boxes, -1);
            LocalizationGraph.LocalizeOGlycan(localizationGraph, scans.First(), commonParameters.ProductMassTolerance, products);
            var pepScore = MetaMorpheusEngine.CalculatePeptideScore(scans.First().TheScan, matchedKnownFragmentIons.Where(p=>p.NeutralTheoreticalProduct.ProductType == ProductType.c || p.NeutralTheoreticalProduct.ProductType == ProductType.zDot).ToList());
            Assert.That((int)pepScore == (int)localizationGraph.TotalScore); //TO THINK: the two score should be the same.

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
            //The data is from Glyco-DIA paper.
            //Get glycanBox
            var glycanBox = OGlycanBoxes[0];

            //Get unmodified peptide, products, allPossible modPos and all boxes.
            Protein protein = new Protein("AATVGSLAGQPLQER", "P16150");
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();
            List<Product> products = new List<Product>();
            peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);

            int[] modPos = GlycoPeptides.GetPossibleModSites(peptide, new string[] { "S", "T" }).OrderBy(p => p).ToArray();
            var boxes = GlycanBox.BuildChildGlycanBoxes( glycanBox.ModIds, GlycanBox.GlobalOGlycans, GlycanBox.GlobalOGlycanMods).ToArray();

            //Load scan.
            CommonParameters commonParameters = new CommonParameters(dissociationType: DissociationType.ETD, trimMsMsPeaks: false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\181217_Fusion_(LC2)_NewObj_Serum_deSA_Jacalin_HRM_4h_ETD_HCD_DDA_mz(400_1200)_21707.mgf");
            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
            var scans = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).ToArray();

            //Known peptideWithMod match.
            var peptideWithMod = GlycoPeptides.GlyGetTheoreticalPeptide(new int[1] {4}, peptide, glycanBox, GlycanBox.GlobalOGlycanMods);
            Assert.That(peptideWithMod.FullSequence == "AAT[O-Glycosylation:N1 on X]VGSLAGQPLQER");
            //List<Product> knownProducts = peptideWithMod.Fragment(DissociationType.EThcD, FragmentationTerminus.Both).ToList();
            List<Product> knownProducts = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.ETD, peptide, peptideWithMod);
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

            //Graph Localization Mod
            string[] modMotifs = new string[modPos.Length];
            int ip = 0;
            foreach (var n in modPos)
            {
                modPos[ip] = n;
                modMotifs[ip] = "S/T";
                ip++;
            }
            LocalizationGraph localizationGraphMod = new LocalizationGraph(modPos, modMotifs, glycanBox, boxes);

            LocalizationGraph.LocalizeMod(localizationGraphMod, scans.First(), commonParameters.ProductMassTolerance, products, GlycoPeptides.GetLocalFragmentGlycan, GlycoPeptides.GetUnlocalFragmentGlycan);
        }

        [Test]
        public static void OGlycoTest_Run()
        {
            var task = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/GlycoSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);

            Directory.CreateDirectory(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"));
            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/P16150.fasta"), false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_4565.mgf");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) }, new List<string> { spectraFile }, new List<DbForTask> { db }, Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData")).Run();
            Directory.Delete(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"), true);
        }

        [Test]
        public static void OGlycoTest_Run2()
        {
            GlycoSearchTask task = new GlycoSearchTask()
            {
                CommonParameters = new CommonParameters
                (
                    dissociationType: DissociationType.ETD
                )
            };
            task._glycoSearchParameters.OxoniumIonFilt = false;

            Directory.CreateDirectory(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"));
            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/P02649.fasta"), false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\181217_Fusion_(LC2)_NewObj_Serum_deSA_Jacalin_HRM_4h_ETD_HCD_DDA_mz(400_1200)_21707.mgf");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) }, new List<string> { spectraFile }, new List<DbForTask> { db }, Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData")).Run();
            Directory.Delete(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"), true);
        }


        [Test]
        public static void OGlycoTest_LocalizeMod()
        {
            //This test proves that the LocalizeMod method is the same as LocalizeOGlycan. 

            //Get glycanBox
            var glycanBox = OGlycanBoxes[19];

            //Get unmodified peptide, products, allPossible modPos and all boxes.
            Protein protein = new Protein("TTGSLEPSSGASGPQVSSVK", "P16150");
            var peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();
            List<Product> products = new List<Product>();
            peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);

            int[] modPos = GlycoPeptides.GetPossibleModSites(peptide, new string[] { "S", "T" }).OrderBy(v => v).ToArray();
            var boxes = GlycanBox.BuildChildGlycanBoxes(glycanBox.ModIds, GlycanBox.GlobalOGlycans, GlycanBox.GlobalOGlycanMods).ToArray();
            Assert.That(boxes.Count() == 6);

            //Get Unlocal Fragment
            var unlocalCost = GlycoPeptides.GetUnlocalFragment(products, modPos, glycanBox);
            Assert.That(unlocalCost.Count == 5); //Basicly, the unlocal are c/z ions that don't localize glycosylation. 

            //Get scan
            CommonParameters commonParameters = new CommonParameters(dissociationType: DissociationType.EThcD, trimMsMsPeaks: false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_4565.mgf");
            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
            var scans = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).ToArray();

            //Known peptideWithMod match.
            var peptideWithMod = GlycoPeptides.GlyGetTheoreticalPeptide(new int[3] { 10, 2, 3 }, peptide, glycanBox, GlycanBox.GlobalOGlycanMods);
            Assert.That(peptideWithMod.FullSequence == "T[O-Glycosylation:H1N1 on X]T[O-Glycosylation:H1N1 on X]GSLEPSS[O-Glycosylation:N1 on X]GASGPQVSSVK");
            List<Product> knownProducts = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.EThcD, peptide, peptideWithMod);
            var matchedKnownFragmentIons = MetaMorpheusEngine.MatchFragmentIons(scans.First(), knownProducts, commonParameters);

            //Graph Localization
            LocalizationGraph localizationGraph = new LocalizationGraph(modPos, glycanBox, boxes, -1);
            LocalizationGraph.LocalizeOGlycan(localizationGraph, scans.First(), commonParameters.ProductMassTolerance, products);
            var allPaths = LocalizationGraph.GetAllHighestScorePaths(localizationGraph.array, localizationGraph.ChildModBoxes);
            var knowPath = new int[8] { 2, 4, 4, 4, 5, 5, 5, 5 };
            Assert.That(Enumerable.SequenceEqual(knowPath, allPaths[0]));

            //LocalizeMod
            var boxMotifs = new string[] { "S/T","S/T", "S/T", "S/T", "S/T", "S/T", "S/T", "S/T" };
            LocalizationGraph localizationGraph0 = new LocalizationGraph(modPos, boxMotifs, glycanBox, boxes, -1);
            //LocalizationGraph.LocalizeMod(localizationGraph0, scans.First(), commonParameters.ProductMassTolerance, 
            //    products.Where(v => v.ProductType == ProductType.c || v.ProductType == ProductType.zDot).ToList(),
            //    GlycoPeptides.GetLocalFragment, GlycoPeptides.GetUnlocalFragment);
            LocalizationGraph.LocalizeMod(localizationGraph0, scans.First(), commonParameters.ProductMassTolerance,
    products.Where(v => v.ProductType == ProductType.c || v.ProductType == ProductType.zDot).ToList(),
    GlycoPeptides.GetLocalFragmentGlycan, GlycoPeptides.GetUnlocalFragmentGlycan);
            var allPaths0 = LocalizationGraph.GetAllHighestScorePaths(localizationGraph0.array, localizationGraph0.ChildModBoxes);
            var knowPath0 = new int[8] { 2, 4, 4, 4, 5, 5, 5, 5 };
            Assert.That(Enumerable.SequenceEqual(knowPath0, allPaths0[0]));

            Assert.That(localizationGraph.TotalScore < localizationGraph0.TotalScore + 0.00001 && localizationGraph.TotalScore > localizationGraph0.TotalScore - 0.00001);

        }

    }
}
