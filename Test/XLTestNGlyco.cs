using Chemistry;
using EngineLayer;
using EngineLayer.CrosslinkSearch;
using EngineLayer.GlycoSearch;
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
using MathNet.Numerics.LinearRegression;
using MathNet.Numerics.LinearAlgebra;

namespace Test
{
    [TestFixture]
    public class XLTestNGlyco
    {

        private static TestDataFile myMsDataFile { get; set; }

        [OneTimeSetUp]
        public static void Setup()
        {
            myMsDataFile = new TestDataFile();
        }


        [Test]
        public static void GlyTest_ModificationSites()
        {
            PeptideWithSetModifications pep = new PeptideWithSetModifications("ELNPTPNVEVNVECR", null); 
            string[] motifs = new string[] { "Nxs", "Nxt"};
            var sites = GlycoPeptides.GetPossibleModSites(pep, motifs);
            Assert.That(sites.Count() == 1 && sites[0] == 4);

            ModificationMotif.TryGetMotif("C", out ModificationMotif motif1);
            Modification mod1 = new Modification(_originalId: "Carbamidomethyl of C", _modificationType: "Common Fixed", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 57.02146372068994);
            ModificationMotif.TryGetMotif("N", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Test of N", _modificationType: "Common Fixed", _target: motif2, _locationRestriction: "Anywhere.");
            var testN = new PeptideWithSetModifications("C[Common Fixed:Carbamidomethyl of C]N[Common Fixed:Test of N]SSDQPKL[Common Fixed:Carbamidomethyl of C]NLSGIETP", new Dictionary<string, Modification> { { "Carbamidomethyl of C", mod1 }, { "Test of N", mod2 } });
            var testSites = GlycoPeptides.GetPossibleModSites(testN, motifs);
            Assert.That(testSites.Count() == 1 && testSites[0] == 11);


            var testC = new PeptideWithSetModifications("TELAAYLSC[Common Fixed:Carbamidomethyl of C]NATK", new Dictionary<string, Modification> { { "Carbamidomethyl of C", mod1 }});
            var testCSites = GlycoPeptides.GetPossibleModSites(testC, motifs);
            Assert.That(testCSites.Count() == 1 && testSites[0] == 11);
        }

        [Test]
        public static void GlyTest_GlyGetTheoreticalFragments()
        {
            //Create Target peptide
            Protein pep = new Protein("TKPREEQYNSTYR", "accession");
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 7);
            var aPeptideWithSetModifications = pep.Digest(digestionParams, new List<Modification>(), new List<Modification>());
            var peptide = aPeptideWithSetModifications.Last();

            //Create glycans
            string[] motifs = new string[] { "Nxs", "Nxt" };
            var sites = GlycoPeptides.GetPossibleModSites(peptide, motifs);
            Glycan glycan = Glycan.Struct2Glycan("(N(F)(N(H(H(N))(H(N)))))", 0);
            var glycanMod = Glycan.NGlycanToModification(glycan);
            Glycan[] glycans = new Glycan[1] { glycan };
            Modification[] modifications = new Modification[1] { glycanMod };

            //Create Scan.
            CommonParameters commonParameters = new CommonParameters(deconvolutionMassTolerance: new PpmTolerance(20), trimMsMsPeaks: false);
            ////The data GlycoTestData/Glyco_3383.mgf is sliced from antibody MSQC4 N-glycopeptide study: 02-06-17_MSQC4 In-Solution Trypsin Digestion
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/Glyco_3383.mgf"); //"25170.mgf"
            MyFileManager myFileManager = new MyFileManager(true);
            var msDataFile = myFileManager.LoadFile(filePath, commonParameters);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(msDataFile, filePath, commonParameters).ToArray();

            //Generate theoretical glyco peptide and theoretical fragment ions.
            var glycopep = GlycoPeptides.GlyGetTheoreticalPeptide(sites.ToArray(), peptide, new Modification[]{ glycanMod });
            List<Product> fragmentIons = GlycoPeptides.GlyGetTheoreticalFragments(GlycoType.NGlycoPep, DissociationType.HCD, peptide, glycopep, sites, glycans);
            List<Product> etdFragmentIons = GlycoPeptides.GlyGetTheoreticalFragments(GlycoType.NGlycoPep, DissociationType.ETD, peptide, glycopep, sites, glycans);
            List<Product> ethcdFragmentIons = GlycoPeptides.GlyGetTheoreticalFragments(GlycoType.NGlycoPep, DissociationType.EThcD, peptide, glycopep, sites, glycans);
            Assert.That(fragmentIons.Count == 63);
            Assert.That(etdFragmentIons.Count == 36);
            Assert.That(ethcdFragmentIons.Count == 87);

            var glycanYIons0 = GlycoPeptides.GetGlycanYIons(aPeptideWithSetModifications.Last(), glycans);
            Assert.That(glycanYIons0.Count == 17);
            var matchedGlycanYIons = MetaMorpheusEngine.MatchFragmentIons(listOfSortedms2Scans[0], glycanYIons0, commonParameters);
            Assert.AreEqual(matchedGlycanYIons.Count, 14);
            var matchedGlycanYIons_MultiCharge = GlycoPeptides.GlyMatchOriginFragmentIons(listOfSortedms2Scans[0], glycanYIons0, commonParameters);
            Assert.AreEqual(matchedGlycanYIons_MultiCharge.Count, 17);

            var YIonIndicator = GlycoPeptides.GetIndicatorYIon(peptide.MonoisotopicMass, "NN");
            var match = GlycoPeptides.MatchIndicatorYIon(listOfSortedms2Scans[0], YIonIndicator, commonParameters);
            Assert.That(match == true);

            //Match Fragment ions.
            var matchedFragmentIons = MetaMorpheusEngine.MatchFragmentIons(listOfSortedms2Scans[0], fragmentIons, commonParameters);
            Assert.That(matchedFragmentIons.Count == 24);

            var coreIons = GlycoPeptides.ScanGetTrimannosylCore(matchedFragmentIons, GlycoType.NGlycoPep);
            Assert.AreEqual(coreIons.Count, 7); 
            var filter = GlycoPeptides.ScanTrimannosylCoreFilter(matchedFragmentIons, GlycoType.NGlycoPep);
            Assert.AreEqual(filter, true);

            //Test new score functions.
            var pg3score = GlycoPeptides.CalculateGlycoPeptideScore(matchedFragmentIons, fragmentIons, commonParameters);

            //Use Graph Localization method
            int[] n_modPos = GlycoPeptides.GetPossibleModSites(peptide, new string[] { "Nxt", "Nxs" }).OrderBy(p => p).ToArray();
            int[] modPos = new int[n_modPos.Length];
            string[] modMotifs = new string[n_modPos.Length];
            int ip = 0;
            foreach (var n in n_modPos)
            {
                modPos[ip] = n;
                modMotifs[ip] = "Nxs/t";
                ip++;
            }
            Array.Sort(modPos, modMotifs);

            //Build GlycanBox
            int[] ids = new int[1] { 0 };
            GlycanBox glycanBox = new GlycanBox(ids, modMotifs, glycans, modifications);
            glycanBox.ChildGlycanBoxes = GlycanBox.BuildChildGlycanBoxes(glycanBox.ModIds, glycans, modifications).ToArray();

            //run localizationGraph
            LocalizationGraph localizationGraph = new LocalizationGraph(modPos, modMotifs, glycanBox, glycanBox.ChildGlycanBoxes);
            List<Product> hcdProducts = new List<Product>();
            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, hcdProducts);

            LocalizationGraph.LocalizeMod(localizationGraph, listOfSortedms2Scans[0], commonParameters.ProductMassTolerance, hcdProducts,
                GlycoPeptides.GetLocalFragmentGlycan, GlycoPeptides.GetUnlocalFragmentGlycan);
        }

        [Test]
        public static void GlyTest_RunTask()
        {
            GlycoSearchTask task = new GlycoSearchTask();
            task._glycoSearchParameters.GlycoSearchType = GlycoSearchType.NGlycanSearch;
            Directory.CreateDirectory(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"));
            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/Q9C0Y4.fasta"), false);
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/yeast_glycan_25170.mgf");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) }, new List<string> { raw }, new List<DbForTask> { db }, Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData")).Run();
            Directory.Delete(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"), true);
        }

        [Test]
        public static void GlyTest_AIETD()
        {
            Protein pep = new Protein("TNSSFIQGFVDHVKEDCDR", "accession");
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 19);
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Carbamidomethyl of C", _modificationType: "Common Fixed", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 57.02146372068994);
            var fixedModifications = new List<Modification>() { mod2 };
            var aPeptideWithSetModifications = pep.Digest(digestionParams, fixedModifications, new List<Modification>());

            string[] motifs = new string[] { "Nxs", "Nxt" };
            var sites = GlycoPeptides.GetPossibleModSites(aPeptideWithSetModifications.Last(), motifs);
            Glycan glycan = Glycan.Struct2Glycan("(N(N(H(H(H(H)))(H(H(H(H(H))))))))", 0);

            Tolerance tolerance = new PpmTolerance(20);
            CommonParameters commonParameters = new CommonParameters(doPrecursorDeconvolution:false, trimMsMsPeaks:false, dissociationType:DissociationType.EThcD, productMassTolerance: tolerance);
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/11901_AIETD.mgf"); 
            MyFileManager myFileManager = new MyFileManager(true);
            var msDataFile = myFileManager.LoadFile(filePath, commonParameters);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(msDataFile, filePath, commonParameters).ToArray();

            var XLPrecusorSearchMode = new SinglePpmAroundZeroSearchMode(commonParameters.PrecursorMassTolerance.Value);
            var precusorMatched = XLPrecusorSearchMode.Accepts(aPeptideWithSetModifications.Last().MonoisotopicMass + (double)glycan.Mass/1E5, listOfSortedms2Scans[0].PrecursorMass);
            Assert.AreEqual(precusorMatched, 0);

            var glycanMod = Glycan.NGlycanToModification(glycan);
            var glycopep = GlycoPeptides.GlyGetTheoreticalPeptide(sites.ToArray(), aPeptideWithSetModifications.Last(), new Modification[] { glycanMod });
            List<Product> fragmentIons = new List<Product>();
            glycopep.Fragment(DissociationType.EThcD, FragmentationTerminus.Both, fragmentIons);
               
            var matchedFragmentIons = MetaMorpheusEngine.MatchFragmentIons(listOfSortedms2Scans[0], fragmentIons, commonParameters);

            using (StreamWriter output = new StreamWriter(Path.Combine(TestContext.CurrentContext.TestDirectory, "11091_NGlyco_AIETD.tsv")))
            {
                foreach (var product in fragmentIons)
                {
                    output.WriteLine(product.Annotation + "\t" + ((double)glycan.Mass / 1E5 - product.NeutralLoss).ToString() + "\t" + product.NeutralMass.ToString());
                }
            }

        }

        [Test]
        public static void GlyTest_OxoniumIons()
        {
            CommonParameters commonParameters = new CommonParameters();
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/yeast_glycan_25170.mgf");
            MyFileManager myFileManager = new MyFileManager(true);
            var msDataFile = myFileManager.LoadFile(filePath, commonParameters);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(msDataFile, filePath, commonParameters).ToArray();
            //Tips: Using debug mode to check the number of oxoniumIons, in this case will be 9.
            MassDiffAcceptor massDiffAcceptor = new SinglePpmAroundZeroSearchMode(20);
            var oxoinumIonsExist = GlycoPeptides.ScanOxoniumIonFilter(listOfSortedms2Scans[0], massDiffAcceptor, commonParameters.DissociationType);
            Assert.AreEqual(oxoinumIonsExist.Where(p=>p>0).Count(), 9);
        }

        [Test]
        public static void GlyTest_BinarySearch()
        {
            //This is just to test how binary search works.
            double[] array = new double[] { 3.44, 3.45, 4.55, 4.55, 4.55, 4.55, 4.55, 4.55, 4.55, 5.66 };
            double x = 3.43;
            double x1 = 3.44;
            double x2 = 3.441;
            double x3 = 3.45;
            double y = 4.44;
            double z = 5.67;
            double d = 4.55;
            double t = 4.56;
            double t1 = 5.66;
            var xid = GlycoPeptides.BinarySearchGetIndex(array, x);
            var xid1 = GlycoPeptides.BinarySearchGetIndex(array, x1);
            var xid2 = GlycoPeptides.BinarySearchGetIndex(array, x2);
            var xid3 = GlycoPeptides.BinarySearchGetIndex(array, x3);

            var yid = GlycoPeptides.BinarySearchGetIndex(array, y);
            var zid = GlycoPeptides.BinarySearchGetIndex(array, z);
            var did = GlycoPeptides.BinarySearchGetIndex(array, d);
            var tid = GlycoPeptides.BinarySearchGetIndex(array, t);
            var tid1 = GlycoPeptides.BinarySearchGetIndex(array, t1);

            Assert.AreEqual(xid, 0);
            Assert.AreEqual(yid, 2);
            Assert.AreEqual(zid, 10); //Index out range
            Assert.AreEqual(did, 2);
            Assert.AreEqual(tid, 9);          
        }         

        [Test]
        public static void GlyTest_WeightedRegression()
        {
            var x = new[] { new[] { 1.0, 4.0 }, new[] { 2.0, 5.0 }, new[] { 3.0, 2.0 } };
            var y = new[] { 15.0, 20, 10 };
            var w = new[] { new[] { 1, 0, 0 }, new[] { 0, 1, 0 }, new[] { 0, 0, 1 } };
            double[] p = MultipleRegression.QR(x, y, intercept: false);
            Matrix<double> xm = Matrix<double>.Build.Dense(x.Count(), x.First().Length, (i, j) => x[i][j]);
            Matrix<double> wm = Matrix<double>.Build.Dense(x.Count(), w.First().Length, (i, j) => w[i][j]);

            Vector<double> yv = Vector<double>.Build.Dense(y.Count(), (i) => y[i]);
            var pv = WeightedRegression.Weighted(xm, yv, wm);
        }

        [Test]
        public static void GlyTest_GraphLocal()
        {
            //Create Target peptide
            Protein pep = new Protein("DANNTQFQFTSR", "accession");
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 7);
            var aPeptideWithSetModifications = pep.Digest(digestionParams, new List<Modification>(), new List<Modification>());
            var peptide = aPeptideWithSetModifications.Last();

            //Create glycans
            string[] motifs = new string[] { "Nxs", "Nxt" };
            var sites = GlycoPeptides.GetPossibleModSites(peptide, motifs);
            Glycan glycan = Glycan.Struct2Glycan("(N(N(H(H(H))(H(H)))))", 0);
            var glycanMod = Glycan.NGlycanToModification(glycan);
            Glycan[] glycans = new Glycan[1] { glycan };
            Modification[] modifications = new Modification[1] { glycanMod };

            //Create Scan.
            CommonParameters commonParameters = new CommonParameters(deconvolutionMassTolerance: new PpmTolerance(20), trimMsMsPeaks: false);
            ////The data  is sliced from yeast
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/yeast_glycan_25170.mgf"); 
            MyFileManager myFileManager = new MyFileManager(true);
            var msDataFile = myFileManager.LoadFile(filePath, commonParameters);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(msDataFile, filePath, commonParameters).ToArray();

            //Generate theoretical glyco peptide and theoretical fragment ions.
            var glycopep = GlycoPeptides.GlyGetTheoreticalPeptide(sites.ToArray(), peptide, new Modification[] { glycanMod });
            List<Product> fragmentIons = GlycoPeptides.GlyGetTheoreticalFragments(GlycoType.NGlycoPep, DissociationType.HCD, peptide, glycopep, sites, glycans);
            Assert.That(fragmentIons.Count == 51);

            var glycanYIons0 = GlycoPeptides.GetGlycanYIons(aPeptideWithSetModifications.Last(), glycans);
            Assert.That(glycanYIons0.Count == 8);
            var matchedGlycanYIons = MetaMorpheusEngine.MatchFragmentIons(listOfSortedms2Scans[0], glycanYIons0, commonParameters);
            Assert.AreEqual(matchedGlycanYIons.Count, 8);
            var matchedGlycanYIons_MultiCharge = GlycoPeptides.GlyMatchOriginFragmentIons(listOfSortedms2Scans[0], glycanYIons0, commonParameters);
            Assert.AreEqual(matchedGlycanYIons_MultiCharge.Count, 12);

            //Match Fragment ions.
            var matchedFragmentIons = MetaMorpheusEngine.MatchFragmentIons(listOfSortedms2Scans[0], fragmentIons, commonParameters);
            Assert.That(matchedFragmentIons.Count == 34);

            var coreIons = GlycoPeptides.ScanGetTrimannosylCore(matchedFragmentIons, GlycoType.NGlycoPep);
            Assert.AreEqual(coreIons.Count, 7); 
            var filter = GlycoPeptides.ScanTrimannosylCoreFilter(matchedFragmentIons, GlycoType.NGlycoPep);
            Assert.AreEqual(filter, true);

            //Use Graph Localization method
            int[] n_modPos = GlycoPeptides.GetPossibleModSites(peptide, new string[] { "Nxt", "Nxs" }).OrderBy(p => p).ToArray();
            int[] modPos = new int[n_modPos.Length];
            string[] modMotifs = new string[n_modPos.Length];
            int ip = 0;
            foreach (var n in n_modPos)
            {
                modPos[ip] = n;
                modMotifs[ip] = "Nxs/t";
                ip++;
            }
            Array.Sort(modPos, modMotifs);

            //Build GlycanBox
            int[] ids = new int[1] { 0 };
            GlycanBox glycanBox = new GlycanBox(ids, modMotifs, glycans, modifications);
            glycanBox.ChildGlycanBoxes = GlycanBox.BuildChildGlycanBoxes(glycanBox.ModIds, glycans, modifications).ToArray();

            //run localizationGraph
            LocalizationGraph localizationGraph = new LocalizationGraph(modPos, modMotifs, glycanBox, glycanBox.ChildGlycanBoxes);
            List<Product> hcdProducts = new List<Product>();
            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, hcdProducts);

            LocalizationGraph.LocalizeMod(localizationGraph, listOfSortedms2Scans[0], commonParameters.ProductMassTolerance, hcdProducts,
                GlycoPeptides.GetLocalFragmentGlycan, GlycoPeptides.GetUnlocalFragmentGlycan);
            var pepScore = MetaMorpheusEngine.CalculatePeptideScore(listOfSortedms2Scans[0].TheScan, matchedFragmentIons.Where(p=> p.NeutralTheoreticalProduct.ProductType!= ProductType.D && p.NeutralTheoreticalProduct.ProductType != ProductType.M).ToList());
            Assert.That((int)localizationGraph.TotalScore == (int)pepScore); //TO DO: this score should be the same
        }

        internal class TestDataFile : MsDataFile
        {
            //Create fake N-Glyco MS data. Including: Full Scan, HCD-pd-HCD (N-glycopeptide). 
            public TestDataFile() : base(2, new SourceFile(null, null, null, null, null))
            {
                var mz = new double[] { };
                var intensities = new double[] { };
                var ScansHere = new List<MsDataScan>();

                var mz1 = new double[] { 1030.4214, 1030.6719, 1030.9227, 1031.1735, 1031.4242, 1031.6741 };
                var intensities1 = new double[] { 828013.6, 1617673.9, 2045691.8, 1464751.9, 825144.2, 599464.8 };
                var MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
                ScansHere.Add(new MsDataScan(MassSpectrum1, 1, 1, true, Polarity.Positive, 1,
                    new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, null, "scan=1"));

                var MassSpectrum2 = new MzSpectrum(mz, intensities, false);
                ScansHere.Add(new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, 1.0,
                    new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 112, 1.0, null, "scan=2", 1030.421508789063,
                    4, 2045691.75, 1030.421508789063, 2, DissociationType.HCD, null, 1030.421508789063));

                var MassSpectrum3 = new MzSpectrum(mz, intensities, false);
                ScansHere.Add(new MsDataScan(MassSpectrum3, 3, 2, true, Polarity.Positive, 1.0,
                    new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 103, 1.0, null, "scan=3", 1030.421508789063,
                    4, 2045691.75, 1030.421508789063, 2, DissociationType.HCD, 2, 1030.421508789063));

                Scans = ScansHere.ToArray();
            }

            public static string FilePath
            {
                get
                {
                    return "GlycoTestDataFile";
                }
            }

            public static string Name
            {
                get
                {
                    return "GlycoTestDataFile";
                }
            }

            public void ReplaceScanArrays(double[] mz, double[] intensities, int i)
            {
                MzSpectrum massSpectrum = new MzSpectrum(mz, intensities, false);
                Scans[i] = new MsDataScan(massSpectrum, Scans[i].OneBasedScanNumber, Scans[i].MsnOrder, Scans[i].IsCentroid,
                    Scans[i].Polarity, Scans[i].RetentionTime, Scans[i].ScanWindowRange, Scans[i].ScanFilter, Scans[i].MzAnalyzer,
                    massSpectrum.SumOfAllY, Scans[i].InjectionTime, null, Scans[i].NativeId, Scans[i].SelectedIonMZ, Scans[i].SelectedIonChargeStateGuess,
                    Scans[i].SelectedIonIntensity, Scans[i].IsolationMz, Scans[i].IsolationWidth, Scans[i].DissociationType, Scans[i].OneBasedPrecursorScanNumber,
                    Scans[i].SelectedIonMonoisotopicGuessMz);
            }
        }
    }
}
