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
        [Test]
        public static void GlyTest_ModificationSites()
        {
            PeptideWithSetModifications pep = new PeptideWithSetModifications("ELNPTPNVEVNVECR", null); 
            string[] motifs = new string[] { "Nxs", "Nxt"};
            var sites = GlycoSpectralMatch.GetPossibleModSites(pep, motifs);
            Assert.That(sites.Count() == 1 && sites[0] == 4);

            ModificationMotif.TryGetMotif("C", out ModificationMotif motif1);
            Modification mod1 = new Modification(_originalId: "Carbamidomethyl of C", _modificationType: "Common Fixed", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 57.02146372068994);
            ModificationMotif.TryGetMotif("N", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Test of N", _modificationType: "Common Fixed", _target: motif2, _locationRestriction: "Anywhere.");
            var testN = new PeptideWithSetModifications("C[Common Fixed:Carbamidomethyl of C]N[Common Fixed:Test of N]SSDQPKL[Common Fixed:Carbamidomethyl of C]NLSGIETP", new Dictionary<string, Modification> { { "Carbamidomethyl of C", mod1 }, { "Test of N", mod2 } });
            var testSites = GlycoSpectralMatch.GetPossibleModSites(testN, motifs);
            Assert.That(testSites.Count() == 1 && testSites[0] == 11);


            var testC = new PeptideWithSetModifications("TELAAYLSC[Common Fixed:Carbamidomethyl of C]NATK", new Dictionary<string, Modification> { { "Carbamidomethyl of C", mod1 }});
            var testCSites = GlycoSpectralMatch.GetPossibleModSites(testC, motifs);
            Assert.That(testCSites.Count() == 1 && testSites[0] == 11);
        }

        [Test]
        public static void GlyTest_GlyGetTheoreticalFragments()
        {
            Protein pep = new Protein("TKPREEQYNSTYR", "accession");
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 7);
            var aPeptideWithSetModifications = pep.Digest(digestionParams, new List<Modification>(), new List<Modification>());
            var peptide = aPeptideWithSetModifications.Last();
            string[] motifs = new string[] { "Nxs", "Nxt" };
            var sites = GlycoSpectralMatch.GetPossibleModSites(peptide, motifs);
            Glycan glycan = Glycan.Struct2Glycan("(N(F)(N(H(H(N))(H(N)))))", 0);

            
            //using (StreamWriter output = new StreamWriter(Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycanFragmentions.txt")))
            //{
            //    foreach (var product in fragmentIons)
            //    {
            //        foreach (var ion in product.Item2)
            //        {
            //            output.WriteLine(ion.Annotation + "\t" + ion.NeutralLoss.ToString() + "\t" + ion.NeutralMass.ToString());
            //        }
            //    }
            //}

            CommonParameters commonParameters = new CommonParameters(deconvolutionMassTolerance: new PpmTolerance(20), trimMsMsPeaks: false);
            //The data GlycoTestData/Glyco_3383.mgf is sliced from antibody MSQC4 N-glycopeptide study: 02-06-17_MSQC4 In-Solution Trypsin Digestion
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/Glyco_3383.mgf"); //"25170.mgf"
            MyFileManager myFileManager = new MyFileManager(true);
            var msDataFile = myFileManager.LoadFile(filePath, commonParameters);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(msDataFile, filePath, commonParameters).ToArray();

            //var glycanMod = Glycan.NGlycanToModification(glycan);
            var glycopep = GlycoPeptides.GenerateNGlycopeptide(sites[0], peptide, glycan);
            List<Product> fragmentIons = GlycoPeptides.NGlyGetTheoreticalFragments(DissociationType.HCD, peptide, glycopep, glycan);

            List<Product> etdFragmentIons = GlycoPeptides.NGlyGetTheoreticalFragments(DissociationType.ETD, peptide, glycopep, glycan);
            List<Product> ethcdFragmentIons = GlycoPeptides.NGlyGetTheoreticalFragments(DissociationType.EThcD, peptide, glycopep, glycan);
            Assert.That(fragmentIons.Count == 75);
            Assert.That(etdFragmentIons.Count == 36);
            Assert.That(ethcdFragmentIons.Count == 99);

            var glycanYIons0 = GlycoPeptides.GetGlycanYIons(aPeptideWithSetModifications.Last(), glycan);
            Assert.That(glycanYIons0.Count == 17);
            var matchedGlycanYIons = MetaMorpheusEngine.MatchFragmentIons(listOfSortedms2Scans[0], glycanYIons0, commonParameters);
            Assert.AreEqual(matchedGlycanYIons.Count, 14);

            //var glycanYIons = GlycoPeptides.GetGlycanYIons(listOfSortedms2Scans[0].PrecursorMass, glycan);
            //Assert.That(glycanYIons.Count == 17);
            var matchedGlycanYIons_MultiCharge = GlycoPeptides.GlyMatchOriginFragmentIons(listOfSortedms2Scans[0], glycanYIons0, commonParameters);
            Assert.AreEqual(matchedGlycanYIons_MultiCharge.Count, 17);

            //TO DO: The neutroloss is not annotated well.
            var matchedFragmentIons = MetaMorpheusEngine.MatchFragmentIons(listOfSortedms2Scans[0], fragmentIons, commonParameters);
            Assert.That(matchedFragmentIons.Count == 24);

            var coreIons = GlycoPeptides.ScanGetTrimannosylCore(matchedFragmentIons, glycan);
            Assert.AreEqual(coreIons.Count, 7);
            var filter = GlycoPeptides.ScanTrimannosylCoreFilter(matchedFragmentIons, glycan);
            Assert.AreEqual(filter, true);
            //var NGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.NGlycanLocations[0], true, false);
            //var bestGlycans = GlycoPeptides.MatchBestGlycan(listOfSortedms2Scans[0], NGlycans.ToArray(), commonParameters).Where(p => p != null && p.Item2 >= 2).OrderByDescending(p => p.Item2).Take(100).OrderBy(p => p.Item3).ToArray(); ;

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
            var sites = GlycoSpectralMatch.GetPossibleModSites(aPeptideWithSetModifications.Last(), motifs);
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
            var glycopep = GlycoPeptides.GenerateNGlycopeptide(sites[0], aPeptideWithSetModifications.Last(), glycan);
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
    }
}
