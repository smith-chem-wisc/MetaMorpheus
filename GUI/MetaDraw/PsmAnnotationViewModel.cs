using Chemistry;
using EngineLayer;
using iTextSharp.text.pdf;
using MassSpectrometry;
using OxyPlot;
using OxyPlot.Annotations;
using OxyPlot.Axes;
using OxyPlot.Series;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Reflection;

namespace ViewModels
{
    public class PsmAnnotationViewModel : INotifyPropertyChanged
    {
        private const double STROKE_THICKNESS_UNANNOTATED = 0.5;
        private const double STROKE_THICKNESS_ANNOTATED = 2.0;
        private PlotModel privateModel;
        private OxyColor variantCrossColor = OxyColors.Green;

        private static Dictionary<ProductType, OxyColor> productTypeDrawColors = new Dictionary<ProductType, OxyColor>
        {
          { ProductType.b, OxyColors.Blue },
          { ProductType.y, OxyColors.Purple },
          { ProductType.c, OxyColors.Gold },
          { ProductType.zDot, OxyColors.Orange },
          { ProductType.D, OxyColors.DodgerBlue },
          { ProductType.M, OxyColors.Firebrick }
        };

        private static Dictionary<ProductType, OxyColor> betaPeptideProductTypeDrawColors = new Dictionary<ProductType, OxyColor>
        {
          { ProductType.b, OxyColors.LightBlue },
          { ProductType.y, OxyColors.MediumPurple },
          { ProductType.c, OxyColors.LightGoldenrodYellow },
          { ProductType.zDot, OxyColors.OrangeRed },
          { ProductType.D, OxyColors.AliceBlue },
          { ProductType.M, OxyColors.LightCoral }
        };

        public PlotModel Model
        {
            get
            {
                return this.privateModel;
            }
            set
            {
                this.privateModel = value;
                NotifyPropertyChanged("Model");
            }
        }

        public event PropertyChangedEventHandler PropertyChanged;

        protected void NotifyPropertyChanged(string propertyName)
        {
            PropertyChangedEventHandler handler = PropertyChanged;
            if (handler != null)
            {
                handler(this, new PropertyChangedEventArgs(propertyName));
            }
        }

        public PsmAnnotationViewModel()
        {
            // Create and Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
            privateModel = new PlotModel { Title = "Spectrum Annotation", Subtitle = "using OxyPlot" };
        }

        // single peptides (not crosslink)
        public void DrawPeptideSpectralMatch(MsDataScan msDataScan, PsmFromTsv psmToDraw,
            bool annotateMz = false, bool annotateCharge = false, int annotateFontSize = 12, bool annotateBold = false)
        {
            // Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
            this.Model = Draw(msDataScan, psmToDraw, annotateMz, annotateCharge, annotateFontSize, annotateBold);
        }

        private PlotModel Draw(MsDataScan msDataScan, PsmFromTsv psmToDraw,
            bool annotateMz = false, bool annotateCharge = false, int annotateFontSize = 12, bool annotateBold = false)
        {
            // x is m/z, y is intensity
            var spectrumMzs = msDataScan.MassSpectrum.XArray;
            var spectrumIntensities = msDataScan.MassSpectrum.YArray;

            string subtitle = psmToDraw.FullSequence;
            if (psmToDraw != null && psmToDraw.BetaPeptideBaseSequence != null)
            {
                subtitle = psmToDraw.FullSequence + "\n" + psmToDraw.BetaPeptideFullSequence;
            }
            PlotModel model = new PlotModel { Title = "Spectrum Annotation of Scan #" + msDataScan.OneBasedScanNumber, DefaultFontSize = 15, Subtitle = subtitle };

            model.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Bottom,
                Title = "m/z",
                Minimum = msDataScan.ScanWindowRange.Minimum,
                Maximum = msDataScan.ScanWindowRange.Maximum,
                AbsoluteMinimum = 0,
                AbsoluteMaximum = msDataScan.ScanWindowRange.Maximum * 2,
                MajorStep = Math.Round(msDataScan.ScanWindowRange.Maximum / 2),
                MinorStep = Math.Round(msDataScan.ScanWindowRange.Maximum / 4)
            });
            model.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Left,
                Title = "Intensity",
                Minimum = 0,
                Maximum = spectrumIntensities.Max(),
                AbsoluteMinimum = 0,
                AbsoluteMaximum = spectrumIntensities.Max() * 2,
                MajorStep = spectrumIntensities.Max(),
                MinorStep = spectrumIntensities.Max() / 2,
                StringFormat = "0e-0",
            });
            model.Axes[1].Zoom(0, spectrumIntensities.Max() * 1.1);

            LineSeries[] allIons = new LineSeries[spectrumMzs.Length];

            // draw the matched peaks; if the PSM is null, we're just drawing the peaks in the scan without annotation, so skip this part
            if (psmToDraw != null)
            {
                List<MatchedFragmentIon> ionsToDraw = new List<MatchedFragmentIon>();

                // check to see if we're drawing a child scan or a parent scan
                if (psmToDraw.Ms2ScanNumber == msDataScan.OneBasedScanNumber)
                {
                    // parent scan
                    ionsToDraw = psmToDraw.MatchedIons;
                }
                else if (msDataScan.OneBasedScanNumber != psmToDraw.Ms2ScanNumber
                    && psmToDraw.ChildScanMatchedIons.Keys.Any(p => p == msDataScan.OneBasedScanNumber))
                {
                    // child scan
                    var scan = psmToDraw.ChildScanMatchedIons.FirstOrDefault(p => p.Key == msDataScan.OneBasedScanNumber);
                    ionsToDraw = scan.Value;
                }

                // annotate peaks (typical PSM, or alpha peptide of CSM)
                foreach (MatchedFragmentIon matchedIon in ionsToDraw)
                {
                    AnnotatePeak(model, allIons, msDataScan, matchedIon, spectrumIntensities, psmToDraw, annotateCharge,
                        annotateMz, annotateBold, annotateFontSize, false);
                }

                // annotate peaks of beta peptide of CSM
                if (psmToDraw.BetaPeptideBaseSequence != null)
                {
                    ionsToDraw = new List<MatchedFragmentIon>();

                    // check to see if we're drawing a child scan or a parent scan
                    if (psmToDraw.Ms2ScanNumber == msDataScan.OneBasedScanNumber)
                    {
                        // parent scan
                        ionsToDraw = psmToDraw.BetaPeptideMatchedIons;
                    }
                    else if (msDataScan.OneBasedScanNumber != psmToDraw.Ms2ScanNumber
                        && psmToDraw.BetaPeptideChildScanMatchedIons.Keys.Any(p => p == msDataScan.OneBasedScanNumber))
                    {
                        // child scan
                        var scan = psmToDraw.BetaPeptideChildScanMatchedIons.FirstOrDefault(p => p.Key == msDataScan.OneBasedScanNumber);
                        ionsToDraw = scan.Value;
                    }

                    foreach (MatchedFragmentIon matchedIon in ionsToDraw)
                    {
                        AnnotatePeak(model, allIons, msDataScan, matchedIon, spectrumIntensities, psmToDraw, annotateCharge,
                            annotateMz, annotateBold, annotateFontSize, true);
                    }
                }
            }

            // draw the remaining unmatched peaks
            for (int i = 0; i < spectrumMzs.Length; i++)
            {
                // peak has already been drawn (it is a matched peak)
                if (allIons[i] != null)
                {
                    continue;
                }

                allIons[i] = new LineSeries();
                allIons[i].Color = OxyColors.DimGray;
                allIons[i].StrokeThickness = STROKE_THICKNESS_UNANNOTATED;
                allIons[i].Points.Add(new DataPoint(spectrumMzs[i], 0));
                allIons[i].Points.Add(new DataPoint(spectrumMzs[i], spectrumIntensities[i]));
                model.Series.Add(allIons[i]);
            }

            // Axes are created automatically if they are not defined      
            return model;
        }

        private PlotModel DrawPdf(MsDataScan msDataScan, PropertyInfo[] properties, PsmFromTsv psm, bool redraw)
        {
            if (redraw)
            {
                this.Model = Draw(msDataScan, psm);
            }

            var y = Model.DefaultYAxis.ActualMaximum - Model.DefaultYAxis.ActualMaximum * 0.03;
            var x = Model.DefaultXAxis.ActualMaximum - Model.DefaultXAxis.ActualMaximum * 0.01;
            var diff = (y - (Model.DefaultYAxis.ActualMaximum * 0.1)) / properties.Length;

            // properties to include
            string[] propertiesToWrite = {
                "Filename",
                "PrecursorCharge",
                "PrecursorMass",
                "PeptideMonoMass",
                "MassDiffDa",
                "MassDiffPpm",
                "Score",
                "DeltaScore",
                "ProteinAccession",
                "ProteinName",
                "GeneName",
                "DecoyContamTarget",
                "QValue",
                "QValueNotch" };

            var propertiesList = properties.Where(p => propertiesToWrite.Contains(p.Name)).OrderBy(p => Array.IndexOf(propertiesToWrite, p.Name)).ToList();

            var displayedProperties = propertiesList.Where(p => p.GetValue(psm) != null); // only display non-null properties

            foreach (PropertyInfo property in displayedProperties)
            {
                // trim property values > 50 characters
                var val = "" + property.GetValue(psm);
                if (val.Length > 50)
                {
                    val = val.Substring(0, 50) + "...";
                }

                var propertyAnnotation = new TextAnnotation
                {
                    Text = property.Name + ": " + val,
                    TextPosition = new DataPoint(x, y),
                    FontSize = 9,
                    StrokeThickness = 0,
                    TextHorizontalAlignment = HorizontalAlignment.Right
                };

                y -= diff;
                Model.Annotations.Add(propertyAnnotation);
            }

            return Model;
        }

        public void DrawPeptideSpectralMatchPdf(MsDataScan msDataScan, PsmFromTsv psm, string fileName, bool redraw)
        {
            var properties = psm.GetType().GetProperties();
            var pdfModel = DrawPdf(msDataScan, properties, psm, redraw);

            string tempPath = Path.Combine(Path.GetDirectoryName(fileName), "sequence.pdf");
            string baseSeqTempPath = Path.Combine(Path.GetDirectoryName(fileName), "annotation.png");

            // exports plot to pdf
            using (var stream = File.Create(tempPath))
            {
                PdfExporter pdf = new PdfExporter { Width = 800, Height = 500 };
                pdf.Export(pdfModel, stream);
            }

            // adds base seq annotation to pdf
            using (Stream inputPdfStream = new FileStream(tempPath, FileMode.Open, FileAccess.Read, FileShare.Read))
            using (Stream inputImageStream = new FileStream(baseSeqTempPath, FileMode.Open, FileAccess.Read, FileShare.Read))
            using (Stream outputPdfStream = new FileStream(fileName, FileMode.Create, FileAccess.Write, FileShare.None))
            {
                var reader = new PdfReader(inputPdfStream);
                var stamper = new PdfStamper(reader, outputPdfStream);
                var pdfContentByte = stamper.GetOverContent(1);

                var image = iTextSharp.text.Image.GetInstance(inputImageStream);
                image.ScaleAbsoluteHeight(500);
                image.ScaleAbsoluteWidth(500);
                image.SetAbsolutePosition(95, 190);
                pdfContentByte.AddImage(image);
                stamper.Close();
            }

            File.Delete(tempPath);
            File.Delete(baseSeqTempPath);
        }

        private void AnnotatePeak(PlotModel model, LineSeries[] allIons, MsDataScan msDataScan, MatchedFragmentIon matchedIon, double[] spectrumIntensities, PsmFromTsv psmToDraw,
            bool annotateCharge, bool annotateMz, bool annotateBold, int annotateFontSize, bool isBetaPeptide)
        {
            OxyColor ionColor;

            if (psmToDraw.VariantCrossingIons.Contains(matchedIon))
            {
                ionColor = variantCrossColor;
            }
            else if (productTypeDrawColors.ContainsKey(matchedIon.NeutralTheoreticalProduct.ProductType))
            {
                ionColor = productTypeDrawColors[matchedIon.NeutralTheoreticalProduct.ProductType];
            }
            else
            {
                ionColor = OxyColors.Turquoise;
            }

            int i = msDataScan.MassSpectrum.GetClosestPeakIndex(matchedIon.NeutralTheoreticalProduct.NeutralMass.ToMz(matchedIon.Charge)).Value;

            // peak line
            allIons[i] = new LineSeries();
            allIons[i].Color = ionColor;
            allIons[i].StrokeThickness = STROKE_THICKNESS_ANNOTATED;
            allIons[i].Points.Add(new DataPoint(matchedIon.Mz, 0));
            allIons[i].Points.Add(new DataPoint(matchedIon.Mz, spectrumIntensities[i]));

            // peak annotation
            string prefix = "";
            if (psmToDraw.BetaPeptideBaseSequence != null)
            {
                if (isBetaPeptide)
                {
                    //prefix = "β";
                    prefix = "B-";
                }
                else
                {
                    //prefix = "α";
                    prefix = "A-";
                }
            }

            string productType = matchedIon.NeutralTheoreticalProduct.ProductType.ToString().ToLower();//.Replace("star", "*").Replace("degree", "°").Replace("dot", "");
            string productNumber = matchedIon.NeutralTheoreticalProduct.TerminusFragment.FragmentNumber.ToString();
            string peakAnnotationText = prefix + productType + productNumber;

            if (matchedIon.NeutralTheoreticalProduct.NeutralLoss != 0)
            {
                peakAnnotationText += "-" + matchedIon.NeutralTheoreticalProduct.NeutralLoss.ToString("F2") + " (" + matchedIon.Mz.ToString("F3") + ")";
            }

            if (annotateCharge)
            {
                peakAnnotationText += "+" + matchedIon.Charge;
            }

            if (annotateMz)
            {
                peakAnnotationText += " (" + matchedIon.Mz.ToString("F3") + ")";
            }

            var peakAnnotation = new TextAnnotation();
            peakAnnotation.Font = "Arial";
            peakAnnotation.FontSize = annotateFontSize;
            peakAnnotation.FontWeight = annotateBold ? FontWeights.Bold : 2.0;
            peakAnnotation.TextColor = ionColor;
            peakAnnotation.StrokeThickness = 0;
            peakAnnotation.Text = peakAnnotationText;
            peakAnnotation.TextPosition = new DataPoint(allIons[i].Points[1].X, allIons[i].Points[1].Y);
            peakAnnotation.TextHorizontalAlignment = HorizontalAlignment.Center;
            model.Annotations.Add(peakAnnotation);

            model.Series.Add(allIons[i]);
        }
    }
}