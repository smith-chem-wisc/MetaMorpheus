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
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using OxyPlot.Annotations;

namespace Test
{
    [TestFixture]
    public static class XLTestNGlyco
    {
        [Test]
        public static void GlyTest_GetKindString()
        {
            string structure = "(N(F)(N(H(H(N))(H(N)))))";
            string kindString = Glycan.GetKindString(structure);
            byte[] kind = new byte[] {3, 4, 0, 0, 1 };
            string kindString2 = Glycan.GetKindString(kind);
            Assert.AreEqual("H3N4A0G0F1", kindString);
            Assert.AreEqual("H3N4A0G0F1", kindString2);
        }

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

            string[] motifs = new string[] { "Nxs", "Nxt" };
            var sites = GlycoSpectralMatch.GetPossibleModSites(aPeptideWithSetModifications.Last(), motifs);
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

            CommonParameters commonParameters = new CommonParameters(deconvolutionMassTolerance: new PpmTolerance(20));
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/Glyco_3383.mgf"); //"25170.mgf"
            MyFileManager myFileManager = new MyFileManager(true);
            var msDataFile = myFileManager.LoadFile(filePath, commonParameters);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(msDataFile, filePath, commonParameters).ToArray();

            var glycanMod = Glycan.NGlycanToModification(glycan);
            var glycopep = GlycoPeptides.GenerateGlycopeptide(sites[0], aPeptideWithSetModifications.Last(), glycan);
            var fragmentIons = glycopep.Fragment(DissociationType.HCD, FragmentationTerminus.Both).ToList();

            var glycanYIons = GlycoPeptides.GetGlycanYIons(listOfSortedms2Scans[0].PrecursorMass, glycan);
            var matchedGlycanYIons = MetaMorpheusEngine.MatchFragmentIons(listOfSortedms2Scans[0], glycanYIons, commonParameters);
            Assert.AreEqual(matchedGlycanYIons.Count, 16);

            var matchedGlycanYIons2 = MetaMorpheusEngine.MatchOriginFragmentIons(listOfSortedms2Scans[0], glycanYIons, commonParameters);
            Assert.AreEqual(matchedGlycanYIons.Count, 17);
            //TO DO: The neutroloss is not annotated well.
            var matchedFragmentIons = MetaMorpheusEngine.MatchFragmentIons(listOfSortedms2Scans[0], fragmentIons, commonParameters);

            var coreIons = GlycoPeptides.ScanGetTrimannosylCore(matchedFragmentIons, glycan);
            Assert.AreEqual(coreIons.Count, 9);
            var filter = GlycoPeptides.ScanTrimannosylCoreFilter(matchedFragmentIons, glycan);
            Assert.AreEqual(filter, true);
            var NGlycans = Glycan.LoadGlycan(GlobalVariables.NGlycanLocation);
            var bestGlycans = GlycoPeptides.MatchBestGlycan(listOfSortedms2Scans[0], NGlycans.ToArray(), commonParameters).Where(p => p != null && p.Item2 >= 2).OrderByDescending(p => p.Item2).Take(100).OrderBy(p => p.Item3).ToArray(); ;

            //Please keep the draw functions, they are important to debug visually.
            DrawPeptideSpectralMatch(listOfSortedms2Scans[0].TheScan, matchedGlycanYIons2, pep.BaseSequence);
            //DrawPeptideSpectralMatch(listOfSortedms2Scans[0].TheScan, matchedGlycanYIons, pep.BaseSequence);
        }

        [Test]
        public static void GlyTest_RunTask()
        {
            GlycoSearchTask task = new GlycoSearchTask();
            task._glycoSearchParameters.IsOGlycoSearch = false;
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
            var glycopep = GlycoPeptides.GenerateGlycopeptide(sites[0], aPeptideWithSetModifications.Last(), glycan);
            var fragmentIons = glycopep.Fragment(DissociationType.EThcD, FragmentationTerminus.Both).ToList();
               
            var matchedFragmentIons = MetaMorpheusEngine.MatchFragmentIons(listOfSortedms2Scans[0], fragmentIons, commonParameters);

            using (StreamWriter output = new StreamWriter(Path.Combine(TestContext.CurrentContext.TestDirectory, "11091_NGlyco_AIETD.tsv")))
            {
                foreach (var product in fragmentIons)
                {
                    output.WriteLine(product.Annotation + "\t" + ((double)glycan.Mass / 1E5 - product.NeutralLoss).ToString() + "\t" + product.NeutralMass.ToString());
                }
            }

            DrawPeptideSpectralMatch(listOfSortedms2Scans[0].TheScan, matchedFragmentIons, pep.BaseSequence);

        }

        [Test]
        public static void GlyTest_OxoniumIons()
        {
            CommonParameters commonParameters = new CommonParameters();
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/yeast_glycan_25170.mgf");
            MyFileManager myFileManager = new MyFileManager(true);
            var msDataFile = myFileManager.LoadFile(filePath, commonParameters);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(msDataFile, filePath, commonParameters).ToArray();
            //Tips: Using debug mode to check the number of oxoniumIons, in this case will be 7.
            MassDiffAcceptor massDiffAcceptor = new SinglePpmAroundZeroSearchMode(20);
            var oxoinumIonsExist = GlycoPeptides.ScanOxoniumIonFilter(listOfSortedms2Scans[0], massDiffAcceptor, commonParameters.DissociationType);
            Assert.AreEqual(5, oxoinumIonsExist);
        }

        [Test]
        public static void GlyTest_DistinguishGlycans()
        {
            Glycan glycan = Glycan.Struct2Glycan("(N(N(H(H(H(H)))(H(H(H(H)))(H(H(H)))))))", 0);
            Glycan glycan2 = Glycan.Struct2Glycan("(N(N(H(H(H))(H(H(H))(H(H(H(H(H)))))))))", 0);

            var test = Glycan.DistingushGlycans(glycan, glycan2);
            Assert.AreEqual(test, true);

            //TO DO: Test the glycan ions. 
            Glycan glycan3 = Glycan.Struct2Glycan("(N(F)(N(H(H(N(H(N(H(N(H))))))(N(H(N(H(N(F)(H(G))))))))(H(N(H(N(H(N(H(A)))))))(N(F)(H(N(F)(H(N(H)(F))))))))))", 8086);
            Glycan glycan4 = Glycan.Struct2Glycan("(N(F)(N(H(H(N(H(N(H(N(H))))))(N(H(N(H(N(F)(H(A))))))))(H(N(H(N(H(N(H(G)))))))(N(F)(H(N(F)(H(N(H)(F))))))))))", 8087);

        }

        [Test]
        public static void GlyTest_BisectHexNAc()
        {
            //The node here is for check the structure of the glycan. 
            Node node = Glycan.Struct2Node("(N(N(H(N)(H(N)(N))(H(N(H))))))"); //This glycan has a bisect hexnac 
            Glycan glycan = Glycan.Struct2Glycan("(N(N(H(N)(H(N)(N))(H(N(H))))))", 0);
            Assert.AreEqual(glycan.Ions.Count, 17);
        }

        [Test]
        public static void GlyTest_GlycanDecoy()
        {
            Glycan glycan = Glycan.Struct2Glycan("(N(N(H(N)(H(N)(N))(H(N(H))))))", 0);
            var test = Glycan.BuildTargetDecoyGlycans(new Glycan[] { glycan});
            Assert.AreEqual(test.Last().Decoy, true);
            foreach (var ion in test.Last().Ions)
            {
                Assert.AreEqual(ion.IonMass + ion.LossIonMass, test.Last().Mass);
            }
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
      
        
        //This is not exactly a test. The function is used for N-Glycan database generation. The function maybe useful in the future.
        public static void GlyTest_GenerateDataBase()
        {         
            var NGlycans = Glycan.LoadGlycan(GlobalVariables.NGlycanLocation);
            string aietdpath = @"GlycoTestData/ComboGlycanDatabase.csv";
            var glycans = Glycan.LoadKindGlycan(aietdpath, NGlycans).ToList();

            string aietdpathWritePath = "E:\\MassData\\Glycan\\GlycanDatabase\\AIETD\\GlycansAIETD.tsv";
            using (StreamWriter output = new StreamWriter(aietdpathWritePath))
            {
                foreach (var glycan in glycans)
                {
                    output.WriteLine(glycan.Mass.ToString() + "\t" + glycan.Struc + "\t" + Glycan.GetKindString(glycan.Kind) + "\t" + Glycan.GetKindString(Glycan.GetKind(glycan.Struc)));
                }
            }
        }

        //This is not exactly a test. The function is used for N-Glycan database generation. The function maybe useful in the future.
        [Test]
        public static void GlyTest_GenerateUniprotDataBase()
        {
            var NGlycans = Glycan.LoadGlycan(GlobalVariables.NGlycanLocation).GroupBy(p => Glycan.GetKindString(p.Kind)).ToDictionary(p => p.Key, p => p.ToList());
            string pathWritePath = "E:\\MassData\\Glycan\\GlycanDatabase\\UniprotGlycanDatabase.txt";


            using (StreamWriter output = new StreamWriter(pathWritePath))
            {
                foreach (var glycan in NGlycans)
                {
                    var mod = Glycan.NGlycanToModification(glycan.Value.First());
                    List<string> temp = new List<string>();
                    temp.Add(mod.ToString());
                    temp.Add(@"//");
                    foreach (var v in temp)
                    {
                        output.WriteLine(v);
                    }                  
                }
            }

        }

        [Test]
        public static void GlyTest_GetAllIonMassFromKind()
        {

            var NGlycans = Glycan.LoadGlycan(GlobalVariables.NGlycanLocation);
            var groupedGlycans = NGlycans.GroupBy(p => Glycan.GetKindString(p.Kind)).ToDictionary(p => p.Key, p => p.ToList());

            byte[] k36101 = new byte[] { 3, 6, 1, 0, 1 };
            var glycan36101 = Glycan.GetAllIonMassFromKind(k36101, groupedGlycans);

            Glycan glycan = Glycan.Struct2Glycan("(N(F)(N(H(H(N))(H(N)))))", 0);
            var x = Glycan.GetAllIonMassFromKind(glycan.Kind, groupedGlycans);
            Assert.AreEqual(x.First().Ions.Count, 14);

            string aietdpath = @"GlycoTestData/ComboGlycanDatabase.csv";
            var glycans = Glycan.LoadKindGlycan(aietdpath, NGlycans).ToList();
            Assert.AreEqual(glycans.Count, 182);

        }

        private static Dictionary<ProductType, OxyColor> productTypeDrawColors = new Dictionary<ProductType, OxyColor>
        {
          { ProductType.b, OxyColors.Blue },
          { ProductType.y, OxyColors.Purple },
          { ProductType.c, OxyColors.Gold },
          { ProductType.zDot, OxyColors.Orange },
          { ProductType.D, OxyColors.DodgerBlue },
          { ProductType.M, OxyColors.Firebrick }
        };

        public static void DrawPeptideSpectralMatch(MsDataScan msDataScan, List<MatchedFragmentIon> matchedFragmentIons, String sequence)
        {
            // x is m/z, y is intensity
            var spectrumMzs = msDataScan.MassSpectrum.XArray;
            var spectrumIntensities = msDataScan.MassSpectrum.YArray;

            string subTitle = "---";//sequence;

            PlotModel model = new PlotModel { Title = "Spectrum Annotation of Scan #" + msDataScan.OneBasedScanNumber, DefaultFontSize = 15, Subtitle = subTitle };
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Title = "m/z", Minimum = 0, Maximum = spectrumMzs.Max() * 1.2, AbsoluteMinimum = 0, AbsoluteMaximum = spectrumMzs.Max() * 5 });
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = "Intensity", Minimum = 0, Maximum = spectrumIntensities.Max() * 1.2, AbsoluteMinimum = 0, AbsoluteMaximum = spectrumIntensities.Max() * 1.3 });
            model.Axes[1].Zoom(0, spectrumIntensities.Max() * 1.1);

            LineSeries[] allIons = new LineSeries[spectrumMzs.Length];

            for (int i = 0; i < spectrumMzs.Length; i++)
            {
                allIons[i] = new LineSeries();
                allIons[i].Color = OxyColors.DimGray;
                allIons[i].StrokeThickness = 0.5;
                allIons[i].Points.Add(new DataPoint(spectrumMzs[i], 0));
                allIons[i].Points.Add(new DataPoint(spectrumMzs[i], spectrumIntensities[i]));
                model.Series.Add(allIons[i]);
            }

            foreach (var peak in matchedFragmentIons)
            {
                OxyColor ionColor = productTypeDrawColors[peak.NeutralTheoreticalProduct.ProductType];

                int i = msDataScan.MassSpectrum.GetClosestPeakIndex(peak.Mz).Value;

                // peak line
                allIons[i] = new LineSeries();
                allIons[i].Color = ionColor;
                allIons[i].StrokeThickness = 1;
                allIons[i].Points.Add(new DataPoint(peak.Mz, 0));
                allIons[i].Points.Add(new DataPoint(peak.Mz, peak.Intensity));

                // peak annotation
                var peakAnnotation = new TextAnnotation();
                peakAnnotation.TextRotation = -60;
                peakAnnotation.Font = "Arial";
                peakAnnotation.FontSize = 4;
                peakAnnotation.FontWeight = 1.0;
                peakAnnotation.TextColor = ionColor;
                peakAnnotation.StrokeThickness = 0;
                //string gly = peak.TheoreticalFragmentIon.ProductType == ProductType.None ? "-Hex" : "";
                peakAnnotation.Text = peak.Annotation;
                peakAnnotation.TextPosition = new DataPoint(allIons[i].Points[1].X, allIons[i].Points[1].Y + peakAnnotation.Text.Length * 1.5 / 4);
                peakAnnotation.TextHorizontalAlignment = HorizontalAlignment.Left;
                model.Annotations.Add(peakAnnotation);

                model.Series.Add(allIons[i]);
            }

            // Axes are created automatically if they are not defined

            // Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
            using (var stream = File.Create(Path.Combine(TestContext.CurrentContext.TestDirectory, "annotation.pdf")))
            {
                PdfExporter pdf = new PdfExporter { Width = 500, Height = 210 };
                pdf.Export(model, stream);
            }
        }

    }
}
