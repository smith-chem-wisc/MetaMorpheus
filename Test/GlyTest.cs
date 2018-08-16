using Chemistry;
using EngineLayer;
using EngineLayer.CrosslinkAnalysis;
using EngineLayer.CrosslinkSearch;
using EngineLayer.Indexing;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;
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
        public static void GlyTest_NGlycoSite()
        {
            var commonParameters = new CommonParameters(doPrecursorDeconvolution: false, cIons: true, zDotIons: true, scoreCutoff: 2, digestionParams: new DigestionParams(minPeptideLength: 5));

            //Create databases contain protein.
            var proteinList = new List<Protein>
            { new Protein("DANNTQFQFTSR", "25170"),
                new Protein("DANNSQFQFTSR", "25171"),
                new Protein("DANNCQFQFTSR", "25172"),
                new Protein("DANPTQFQFTSR", "25173"),
                new Protein("DANNCCFQFTSR", "25174"),
                new Protein("DANNCCFQFNSS", "25174"),
                new Protein("DANNCCFQFTNS", "25174")
            };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            ModificationWithMass mod1 = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            ModificationWithMass mod2 = new ModificationWithMass("Carbamidomethyl of C", "Common Fixed", motif2, TerminusLocalization.Any, 57.02146372068994);
            var variableModifications = new List<ModificationWithMass>() { mod1 };
            var fixedModifications = new List<ModificationWithMass>() { mod2 };
            var localizeableModifications = new List<ModificationWithMass>();

            var lp = new List<ProductType> { ProductType.BnoB1ions, ProductType.Y };
            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            foreach (var mod in fixedModifications)
                modsDictionary.Add(mod, 0);
            int i = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }
            foreach (var mod in localizeableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }

            //Generate digested peptide lists.
            List<PeptideWithSetModifications> digestedList = new List<PeptideWithSetModifications>();
            foreach (var item in proteinList)
            {
                var digested = item.Digest(commonParameters.DigestionParams, fixedModifications, variableModifications).ToList();
                digestedList.AddRange(digested);
            }

            foreach (var fdfd in digestedList)
            {
                fdfd.CompactPeptide(TerminusType.None);
            }

            var t1 = PsmCross.NGlyPosCal(digestedList[0].CompactPeptide(TerminusType.None));
            var t2 = PsmCross.NGlyPosCal(digestedList[1].CompactPeptide(TerminusType.None));
            var t3 = PsmCross.NGlyPosCal(digestedList[2].CompactPeptide(TerminusType.None));
            var t4 = PsmCross.NGlyPosCal(digestedList[3].CompactPeptide(TerminusType.None));
            //TO DO: Residue 'C' or 'S', 'T' maybe modified, so that the NGlyPosCal is limited.
            var t5 = PsmCross.NGlyPosCal(digestedList[4].CompactPeptide(TerminusType.None));
            var t6 = PsmCross.NGlyPosCal(digestedList[5].CompactPeptide(TerminusType.None));
            var t7 = PsmCross.NGlyPosCal(digestedList[6].CompactPeptide(TerminusType.None));

            Assert.AreEqual(t1, new List<int>() { 3 });
            Assert.AreEqual(t2, new List<int>() { 3 });
            Assert.AreEqual(t3, new List<int>() { 3 });
            Assert.AreEqual(t4, new List<int>());
            Assert.AreEqual(t5, new List<int>() { 3, 4 });
        }

        [Test]
        public static void GlyTest_OxoniumIons()
        {
            //Get MS2 scans.
            var commonParameters = new CommonParameters(doPrecursorDeconvolution: false, cIons: true, zDotIons: true, scoreCutoff: 2, digestionParams: new DigestionParams(minPeptideLength: 5));
            var myMsDataFile = new GlyTestDataFile();
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, commonParameters.DoPrecursorDeconvolution, commonParameters.UseProvidedPrecursorInfo, commonParameters.DeconvolutionIntensityRatio, commonParameters.DeconvolutionMaxAssumedChargeState, commonParameters.DeconvolutionMassTolerance).ToArray();
            //Tips: Using debug mode to check the number of oxoniumIons, in this case will be 7.
            var oxoinumIonsExist = TwoPassCrosslinkSearchEngine.ScanOxoniumIonFilter(listOfSortedms2Scans[0]);          
            Assert.AreEqual(true, oxoinumIonsExist);
        }

        [Test]
        public static void GlyTest_FragmentIons()
        {
            var commonParameters = new CommonParameters(doPrecursorDeconvolution: false, cIons: true, zDotIons: true, scoreCutoff: 2, digestionParams: new DigestionParams(minPeptideLength: 5));
            var xlSearchParameters = new XlSearchParameters { SearchGlycoWithBgYgIndex = false };

            //Create databases contain protein.
            var proteinList = new List<Protein> { new Protein("DANNTQFQFTSR", "25170") };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            ModificationWithMass mod1 = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            ModificationWithMass mod2 = new ModificationWithMass("Carbamidomethyl of C", "Common Fixed", motif2, TerminusLocalization.Any, 57.02146372068994);
            var variableModifications = new List<ModificationWithMass>() { mod1 };
            var fixedModifications = new List<ModificationWithMass>() { mod2 };
            var localizeableModifications = new List<ModificationWithMass>();

            var lp = new List<ProductType> { ProductType.BnoB1ions, ProductType.Y};
            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            foreach (var mod in fixedModifications)
                modsDictionary.Add(mod, 0);
            int i = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }
            foreach (var mod in localizeableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }

            //Generate digested peptide lists.
            List<PeptideWithSetModifications> digestedList = new List<PeptideWithSetModifications>();
            foreach (var item in proteinList)
            {
                var digested = item.Digest(commonParameters.DigestionParams, fixedModifications, variableModifications).ToList();
                digestedList.AddRange(digested);
            }

            foreach (var fdfd in digestedList)
            {
                fdfd.CompactPeptide(TerminusType.None);
            }

            //Run index engine
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, lp, 1, DecoyType.Reverse, new List<DigestionParams>
            { commonParameters.DigestionParams }, commonParameters, 30000, xlSearchParameters.SearchGlycoWithBgYgIndex, new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();

            var fragmentIndexAll = indexResults.FragmentIndex.Select((s, j) => new { j, s }).Where(p => p.s != null).Select(t => t.j).ToList();
            Assert.IsTrue(fragmentIndexAll.Count() > 0);

            //Get MS2 scans.
            var myMsDataFile = new GlyTestDataFile();
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, commonParameters.DoPrecursorDeconvolution, commonParameters.UseProvidedPrecursorInfo, commonParameters.DeconvolutionIntensityRatio, commonParameters.DeconvolutionMaxAssumedChargeState, commonParameters.DeconvolutionMassTolerance).ToArray();
            
            //TwoPassCrosslinkSearchEngine.Run().
            List<PsmCross> newPsms = new List<PsmCross>();
            new TwoPassCrosslinkSearchEngine(newPsms, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, lp, 0, commonParameters, false, true, false, xlSearchParameters.XlPrecusorMsTl, null, xlSearchParameters.CrosslinkSearchTop, xlSearchParameters.CrosslinkSearchTopNum, xlSearchParameters.XlQuench_H2O, xlSearchParameters.XlQuench_NH2, xlSearchParameters.XlQuench_Tris, xlSearchParameters.XlCharge_2_3, new List<string> { }).Run();

            var compactPeptideToProteinPeptideMatch = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            new CrosslinkAnalysisEngine(newPsms, compactPeptideToProteinPeptideMatch, proteinList, variableModifications, fixedModifications, lp, null, null, TerminusType.None, commonParameters, new List<string> { }).Run();

            DrawPeptideSpectralMatch(myMsDataFile.GetAllScansList()[1], newPsms.First());
        }

        private static Dictionary<ProductType, OxyColor> productTypeDrawColors = new Dictionary<ProductType, OxyColor>
        { { ProductType.B, OxyColors.Blue },
          { ProductType.BnoB1ions, OxyColors.Blue },
          { ProductType.Y, OxyColors.Red },
          { ProductType.C, OxyColors.ForestGreen },
          { ProductType.Zdot, OxyColors.Orange },
         { ProductType.X, OxyColors.Violet },
        { ProductType.None, OxyColors.Purple }};

        public static void DrawPeptideSpectralMatch(MsDataScan msDataScan, PsmCross psmToDraw)
        {
            // x is m/z, y is intensity
            var spectrumMzs = msDataScan.MassSpectrum.XArray;
            var spectrumIntensities = msDataScan.MassSpectrum.YArray;

            string subTitle = psmToDraw.FullSequence;
            if (psmToDraw.BetaPsmCross!= null)
            {
                subTitle += "-" + psmToDraw.BetaPsmCross.FullSequence;
            }
            PlotModel model = new PlotModel { Title = "Spectrum Annotation of Scan #" + msDataScan.OneBasedScanNumber, DefaultFontSize = 15, Subtitle = subTitle };
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Title = "m/z", Minimum = 0, Maximum = spectrumMzs.Max() * 1.2, AbsoluteMinimum = 0, AbsoluteMaximum = spectrumMzs.Max() * 5 });
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = "Intensity", Minimum = 0, Maximum = spectrumIntensities.Max() * 1.2, AbsoluteMinimum = 0, AbsoluteMaximum = spectrumIntensities.Max() * 1.3 });
            model.Axes[1].Zoom(0, spectrumIntensities.Max() * 1.1);

            LineSeries[] allIons = new LineSeries[spectrumMzs.Length];

            // draw the matched peaks; if the PSM is null, we're just drawing the peaks in the scan without annotation, so skip this part
            if (psmToDraw != null)
            {
                if (psmToDraw.BetaPsmCross != null)
                {
                    psmToDraw.MatchedIons = psmToDraw.MatchedIons.Concat(psmToDraw.BetaPsmCross.MatchedIons).ToList();
                }
                foreach (var peak in psmToDraw.MatchedIons)
                {
                    OxyColor ionColor = productTypeDrawColors[peak.TheoreticalFragmentIon.ProductType];

                    int i = msDataScan.MassSpectrum.GetClosestPeakIndex(peak.Mz).Value;

                    // peak line
                    allIons[i] = new LineSeries();
                    allIons[i].Color = ionColor;
                    allIons[i].StrokeThickness = 1;
                    allIons[i].Points.Add(new DataPoint(peak.Mz, 0));
                    allIons[i].Points.Add(new DataPoint(peak.Mz, spectrumIntensities[i]));

                    // peak annotation
                    var peakAnnotation = new TextAnnotation();
                    peakAnnotation.TextRotation = -60;
                    peakAnnotation.Font = "Arial";
                    peakAnnotation.FontSize = 4;
                    peakAnnotation.FontWeight = 1.0;
                    peakAnnotation.TextColor = ionColor;
                    peakAnnotation.StrokeThickness = 0;
                    //string gly = peak.TheoreticalFragmentIon.ProductType == ProductType.None ? "-Hex" : "";
                    peakAnnotation.Text = "(" + peak.Mz.ToString("F3") + "@+"+ peak.TheoreticalFragmentIon.Charge.ToString() + ") " + peak.TheoreticalFragmentIon.ProductType.ToString().ToLower().First() + "-" + peak.TheoreticalFragmentIon.IonNumber;
                    peakAnnotation.TextPosition = new DataPoint(allIons[i].Points[1].X, allIons[i].Points[1].Y + peakAnnotation.Text.Length * 1.5 / 4);
                    peakAnnotation.TextHorizontalAlignment = HorizontalAlignment.Left;
                    model.Annotations.Add(peakAnnotation);

                    model.Series.Add(allIons[i]);
                }           
            }

            for (int i = 0; i < spectrumMzs.Length; i++)
            {
                allIons[i] = new LineSeries();
                allIons[i].Color = OxyColors.DimGray;
                allIons[i].StrokeThickness = 0.5;
                allIons[i].Points.Add(new DataPoint(spectrumMzs[i], 0));
                allIons[i].Points.Add(new DataPoint(spectrumMzs[i], spectrumIntensities[i]));
                model.Series.Add(allIons[i]);
            }

            // Axes are created automatically if they are not defined

            // Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
            using (var stream = File.Create(Path.Combine(TestContext.CurrentContext.TestDirectory, psmToDraw.ScanNumber + ".pdf")))
            {
                PdfExporter pdf = new PdfExporter { Width = 500, Height = 210 };
                pdf.Export(model, stream);
            }
        }
    }

    internal class GlyTestDataFile : MsDataFile
    {
        public GlyTestDataFile() : base(2, new SourceFile(null, null, null, null, null))
        {
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/25170.mgf");
            double[] mz;
            double[] intensity;
            GetScan(raw, out mz, out intensity);

            var mz1 = new double[] { 2644.06992.ToMz(2) };
            var intensities1 = new double[] { 1 };
            var MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
            var ScansHere = new List<MsDataScan> { new MsDataScan(MassSpectrum1, 1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Orbitrap, 1000, 1, null, "scan=1") };

            var MassSpectrum2 = new MzSpectrum(mz, intensity, false);
            ScansHere.Add(new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, 1.0,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Orbitrap, intensity.Sum(), 1.0, null, "scan=2", 2644.06992.ToMz(2),
                2, 1, 2644.06992.ToMz(2), 2, DissociationType.HCD, 1, 2644.06992.ToMz(2)));

            Scans = ScansHere.ToArray();
        }

        public string FilePath
        {
            get
            {
                return "GlyTestDataFile";
            }
        }

        public string Name
        {
            get
            {
                return "GlyTestDataFile";
            }
        }

        public void ReplaceFirstScanArrays(double[] mz, double[] intensities)
        {
            MzSpectrum massSpectrum = new MzSpectrum(mz, intensities, false);
            Scans[0] = new MsDataScan(massSpectrum, Scans[0].OneBasedScanNumber, Scans[0].MsnOrder, Scans[0].IsCentroid, Scans[0].Polarity, Scans[0].RetentionTime, Scans[0].ScanWindowRange, Scans[0].ScanFilter, Scans[0].MzAnalyzer, massSpectrum.SumOfAllY, Scans[0].InjectionTime, null, Scans[0].NativeId);
        }

        private void GetScan(string filePath, out double[] mz, out double[]intensity)
        {
            List<double> mzList = new List<double>();
            List<double> intensityList = new List<double>();
            
            using (StreamReader glycans = new StreamReader(filePath))
            {            
                while (glycans.Peek() != -1)
                {
                    string line = glycans.ReadLine();
                    //TO DO: Read .mgf file in a smarter way.
                    if (!line.StartsWith("B") && !line.StartsWith("T") && !line.StartsWith("R") && !line.StartsWith("P") && !line.StartsWith("T") && !line.StartsWith("C") && !line.StartsWith("E"))
                    {       
                        mzList.Add(Convert.ToDouble(line.Split(null)[0]));
                        intensityList.Add(Convert.ToDouble(line.Split(null)[1]));
                    }
                }
            }
            mz = mzList.ToArray();
            intensity = intensityList.ToArray();
        }
    }
}
