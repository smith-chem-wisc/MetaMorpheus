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

namespace Test
{
    [TestFixture]
    public static class GlyTest
    {
        [Test]
        public static void GlyTest_GetKindString()
        {
            string structure = "(N(F)(N(H(H(N))(H(N)))))";
            string kindString = Glycan.GetKindString(structure);
            int[] kind = new int[] {3, 4, 0, 0, 1 };
            string kindString2 = Glycan.GetKindString(kind);
            Assert.AreEqual("H3N4F1", kindString);
            Assert.AreEqual("H3N4F1", kindString2);
        }

        [Test]
        public static void GlyTest_ModificationSites()
        {
            PeptideWithSetModifications pep = new PeptideWithSetModifications("ELNPTPNVEVNVECR", null); 
            string[] motifs = new string[] { "Nxs", "Nxt"};
            var sites = CrosslinkSpectralMatch.GetPossibleModSites(pep, motifs);
            Assert.That(sites.Count() == 1 && sites[0] == 2);
        }

        [Test]
        public static void GlyTest_GlyGetTheoreticalFragments()
        {
            Protein pep = new Protein("TKPREEQYNSTYR", "accession");
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 7);
            var aPeptideWithSetModifications = pep.Digest(digestionParams, new List<Modification>(), new List<Modification>());

            string[] motifs = new string[] { "Nxs", "Nxt" };
            var sites = CrosslinkSpectralMatch.GetPossibleModSites(aPeptideWithSetModifications.Last(), motifs);
            Glycan glycan = Glycan.Struct2Glycan("(N(F)(N(H(H(N))(H(N)))))", 0);
            var fragmentIons = GlycoPeptides.GlyGetTheoreticalFragments(DissociationType.HCD, sites, aPeptideWithSetModifications.Last(), glycan).ToList();

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
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/Glyco_3383.mgf"); //"XlTestData/25170.mgf"
            MyFileManager myFileManager = new MyFileManager(true);
            var msDataFile = myFileManager.LoadFile(filePath, 300, 0.01, true, true, commonParameters);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(msDataFile, filePath, commonParameters).ToArray();
            //TO DO: The neutroloss is not annotated well.
            var matchedFragmentIons = MetaMorpheusEngine.MatchFragmentIons(listOfSortedms2Scans[0], fragmentIons.Select(p=>p.Item2).ToList().First(), commonParameters);

            DrawPeptideSpectralMatch(listOfSortedms2Scans[0].TheScan, matchedFragmentIons, pep.BaseSequence);
        }

        [Test]
        public static void GlyTest_Topdown()
        {
            Protein pep = new Protein("MAVMAPRTLVLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEPSSQPTIPIVGIIAGLVLFGAVITGAVVAAVMWRRKSSDRKGGSYSQAASSDSAQGSDVSLTACKV", "accession");
            DigestionParams digestionParams = new DigestionParams(protease: "top-down");
            var aPeptideWithSetModifications = pep.Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();

            string[] motifs = new string[] { "Nxs", "Nxt" };
            var sites = CrosslinkSpectralMatch.GetPossibleModSites(aPeptideWithSetModifications, motifs);
            Glycan glycan = Glycan.Struct2Glycan("(N(N(H(H)(H(H)(H)))))", 0);
            var fragmentIons = GlycoPeptides.GlyGetTheoreticalFragments(DissociationType.HCD, sites, aPeptideWithSetModifications, glycan).ToList();

            CommonParameters commonParameters = new CommonParameters(deconvolutionMassTolerance: new PpmTolerance(20));
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/3554.mgf");
            MyFileManager myFileManager = new MyFileManager(true);
            var msDataFile = myFileManager.LoadFile(filePath, 500, 0.01, true, true, commonParameters);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(msDataFile, filePath, commonParameters).ToArray();
            //TO DO: The neutroloss is not annotated well.
            var matchedFragmentIons = MetaMorpheusEngine.MatchFragmentIons(listOfSortedms2Scans[0], fragmentIons.Select(p => p.Item2).ToList().First(), commonParameters);

            DrawPeptideSpectralMatch(listOfSortedms2Scans[0].TheScan, matchedFragmentIons, pep.BaseSequence);
        }
        [Test]
        public static void GlyTest_RunTask()
        {
            var task = Toml.ReadFile<XLSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/GlycoSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            Directory.CreateDirectory(Path.Combine(Environment.CurrentDirectory, @"TESTXlTestData"));
            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/Q9C0Y4.fasta"), false);
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/25170.mgf");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) }, new List<string> { raw }, new List<DbForTask> { db }, Path.Combine(Environment.CurrentDirectory, @"TESTXlTestData")).Run();
            //Directory.Delete(Path.Combine(Environment.CurrentDirectory, @"TESTXlTestData"), true);
        }

        [Test]
        public static void GlyTest_OxoniumIons()
        {
            CommonParameters commonParameters = new CommonParameters();
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/25170.mgf");
            MyFileManager myFileManager = new MyFileManager(true);
            var msDataFile = myFileManager.LoadFile(filePath, 300, 0.01, true, true, commonParameters);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(msDataFile, filePath, commonParameters).ToArray();
            ////Tips: Using debug mode to check the number of oxoniumIons, in this case will be 7.
            var oxoinumIonsExist = CrosslinkSearchEngine.ScanOxoniumIonFilter(listOfSortedms2Scans[0], commonParameters.DissociationType);
            Assert.AreEqual(true, oxoinumIonsExist);
        }

        [Test]
        public static void GlyTest_BinarySearch()
        {
            double[] array = new double[] { 3.44, 3.45, 4.55, 5.66 };
            double x = 3.43;
            double y = 4.44;
            double z = 5.67;
            var xid = Array.BinarySearch(array, x);
            if (xid < 0) { xid = ~xid; }
            var yid = Array.BinarySearch(array, y);
            if (yid < 0) { yid = ~yid; }
            var zid = Array.BinarySearch(array, z);
            if (zid < 0) { zid = ~zid; }
        }

        private static Dictionary<ProductType, OxyColor> productTypeDrawColors = new Dictionary<ProductType, OxyColor>
        {
          { ProductType.b, OxyColors.Blue },
          { ProductType.y, OxyColors.Purple },
          { ProductType.c, OxyColors.Gold },
          { ProductType.zPlusOne, OxyColors.Orange },
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
