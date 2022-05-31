using Chemistry;
using EngineLayer;
using GuiFunctions.MetaDraw;
using iText.IO.Image;
using iText.Kernel.Pdf;
using MassSpectrometry;
using mzPlot;
using OxyPlot;
using OxyPlot.Annotations;
using OxyPlot.Axes;
using OxyPlot.Series;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Media.Imaging;

namespace GuiFunctions
{
    /// <summary>
    /// Class for the peptide spectrum match plot within the metadraw window
    /// </summary>
    public class PeptideSpectrumMatchPlot : Plot
    {
        public MsDataScan Scan { get; protected set; }
        protected PsmFromTsv SpectrumMatch;

        public PeptideSpectrumMatchPlot(OxyPlot.Wpf.PlotView plotView, PsmFromTsv psm, MsDataScan scan,
            List<MatchedFragmentIon> matchedFragmentIons, bool annotateProperties = true, LibrarySpectrum librarySpectrum = null, bool stationarySequence = false) : base(plotView)
        {
            Model.Title = string.Empty;
            Model.Subtitle = string.Empty;
            SpectrumMatch = psm;
            Scan = scan;

            DrawSpectrum();
            AnnotateMatchedIons(isBetaPeptide: false, matchedFragmentIons);

            if (annotateProperties)
            {
                AnnotateProperties();
            }

            ZoomAxes(matchedFragmentIons);

            if (librarySpectrum != null)
            {
                AnnotateLibraryIons(isBetaPeptide: false, librarySpectrum.MatchedFragmentIons);
            }

            RefreshChart();
        }

        public void ExportToPdf(string path, Canvas stationarySequence, double width = 700, double height = 370)
        {
            // exports spectrum annotation w/o base seq annotation
            string tempPdfPath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(path), "temp.pdf");
            string tempPngPath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(path), "annotation.png");
            base.ExportToPdf(tempPdfPath, width, height); // Change path back to temp pdf path

            // scales for desired DPI
            double dpiScale = MetaDrawSettings.CanvasPdfExportDpi / 96.0;

            // save base seq as PNG
            stationarySequence.Measure(new Size((int)stationarySequence.Width + 30, (int)stationarySequence.Height));
            
            stationarySequence.Arrange(new Rect(new Size((int)stationarySequence.Width + 30, (int)stationarySequence.Height)));

            RenderTargetBitmap renderBitmap = new RenderTargetBitmap((int)(dpiScale * stationarySequence.Width + 30), (int)(dpiScale * stationarySequence.Height),
                MetaDrawSettings.CanvasPdfExportDpi, MetaDrawSettings.CanvasPdfExportDpi, PixelFormats.Pbgra32);

            renderBitmap.Render(stationarySequence);

            PngBitmapEncoder encoder = new PngBitmapEncoder();

            encoder.Frames.Add(BitmapFrame.Create(renderBitmap));


            using (FileStream file = File.Create(tempPngPath))
            {
                encoder.Save(file);
            }

            // adds base seq annotation to pdf
            PdfDocument pdfDoc = new PdfDocument(new PdfReader(tempPdfPath), new PdfWriter(path));

            iText.Layout.Document document = new iText.Layout.Document(pdfDoc);

            ImageData imgData = ImageDataFactory.Create(tempPngPath);

            iText.Layout.Element.Image img = new iText.Layout.Element.Image(imgData);
            img.SetMarginLeft((float)(-1.0 * stationarySequence.Margin.Left) + 10);
            img.SetMarginTop(-30);
            img.ScaleToFit((float)stationarySequence.Width, (float)stationarySequence.Height);


            document.Add(img);


            document.Close();
            pdfDoc.Close();

            // delete temp files
            File.Delete(tempPdfPath);
            File.Delete(tempPngPath);
        }

        protected void DrawSpectrum()
        {
            // set up axes
            this.Model.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Bottom,
                Title = "m/z",
                Minimum = Scan.ScanWindowRange.Minimum,
                Maximum = Scan.ScanWindowRange.Maximum,
                AbsoluteMinimum = Math.Max(0, Scan.ScanWindowRange.Minimum - 100),
                AbsoluteMaximum = Scan.ScanWindowRange.Maximum + 100,
                MajorStep = 200,
                MinorStep = 200,
                MajorTickSize = 2,
                TitleFontWeight = OxyPlot.FontWeights.Bold,
                TitleFontSize = 14
            });

            this.Model.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Left,
                Title = "Intensity",
                Minimum = 0,
                Maximum = Scan.MassSpectrum.YArray.Max(),
                AbsoluteMinimum = 0,
                AbsoluteMaximum = Scan.MassSpectrum.YArray.Max() * 2,
                MajorStep = Scan.MassSpectrum.YArray.Max() / 10,
                MinorStep = Scan.MassSpectrum.YArray.Max() / 10,
                StringFormat = "0e-0",
                MajorTickSize = 2,
                TitleFontWeight = OxyPlot.FontWeights.Bold,
                TitleFontSize = 14,
                AxisTitleDistance = 10,
                ExtraGridlines = new double[] { 0 },
                ExtraGridlineColor = OxyColors.Black,
                ExtraGridlineThickness = 1
            });

            // draw all peaks in the scan
            for (int i = 0; i < Scan.MassSpectrum.XArray.Length; i++)
            {
                double mz = Scan.MassSpectrum.XArray[i];
                double intensity = Scan.MassSpectrum.YArray[i];

                DrawPeak(mz, intensity, MetaDrawSettings.StrokeThicknessUnannotated, MetaDrawSettings.UnannotatedPeakColor, null);
            }
        }

        protected void AnnotateMatchedIons(bool isBetaPeptide, List<MatchedFragmentIon> matchedFragmentIons, bool useLiteralPassedValues = false)
        {
            foreach (MatchedFragmentIon matchedIon in matchedFragmentIons)
            {
                AnnotatePeak(matchedIon, isBetaPeptide, useLiteralPassedValues);
            }
        }

        protected void AnnotateLibraryIons(bool isBetaPeptide, List<MatchedFragmentIon> libraryIons)
        {
            // figure out the sum of the intensities of the matched fragment ions
            double sumOfMatchedIonIntensities = 0;
            double sumOfLibraryIntensities = 0;
            foreach (var libraryIon in libraryIons)
            {
                var matchedIon = SpectrumMatch.MatchedIons.FirstOrDefault(p =>
                    p.NeutralTheoreticalProduct.ProductType == libraryIon.NeutralTheoreticalProduct.ProductType
                    && p.NeutralTheoreticalProduct.FragmentNumber == libraryIon.NeutralTheoreticalProduct.FragmentNumber);

                if (matchedIon == null)
                {
                    continue;
                }

                int i = Scan.MassSpectrum.GetClosestPeakIndex(libraryIon.Mz);
                double intensity = Scan.MassSpectrum.YArray[i];
                sumOfMatchedIonIntensities += intensity;
                sumOfLibraryIntensities += libraryIon.Intensity;
            }

            double multiplier = -1 * sumOfMatchedIonIntensities / sumOfLibraryIntensities;

            List<MatchedFragmentIon> mirroredLibraryIons = new List<MatchedFragmentIon>();

            foreach (MatchedFragmentIon libraryIon in libraryIons)
            {
                var neutralProduct = new Product(libraryIon.NeutralTheoreticalProduct.ProductType, libraryIon.NeutralTheoreticalProduct.Terminus,
                    libraryIon.NeutralTheoreticalProduct.NeutralMass, libraryIon.NeutralTheoreticalProduct.FragmentNumber,
                    libraryIon.NeutralTheoreticalProduct.AminoAcidPosition, libraryIon.NeutralTheoreticalProduct.NeutralLoss);

                mirroredLibraryIons.Add(new MatchedFragmentIon(ref neutralProduct, libraryIon.Mz, multiplier * libraryIon.Intensity, libraryIon.Charge));
            }

            AnnotateMatchedIons(isBetaPeptide, mirroredLibraryIons, useLiteralPassedValues: true);

            // zoom to accomodate the mirror plot
            double min = mirroredLibraryIons.Min(p => p.Intensity) * 1.2;
            this.Model.Axes[1].AbsoluteMinimum = min * 2;
            this.Model.Axes[1].AbsoluteMaximum = -min * 2;
            this.Model.Axes[1].Zoom(min, -min);
            this.Model.Axes[1].LabelFormatter = DrawnSequence.YAxisLabelFormatter;
        }

        protected void AnnotatePeak(MatchedFragmentIon matchedIon, bool isBetaPeptide, bool useLiteralPassedValues = false)
        {
            OxyColor ionColor;

            if (SpectrumMatch.VariantCrossingIons.Contains(matchedIon))
            {
                ionColor = MetaDrawSettings.VariantCrossColor;
            }
            else if (matchedIon.NeutralTheoreticalProduct.SecondaryProductType != null) //if internal fragment
            {
                ionColor = OxyColors.Purple;
            }
            else if (isBetaPeptide)
            {
                ionColor = MetaDrawSettings.BetaProductTypeToColor[matchedIon.NeutralTheoreticalProduct.ProductType];
            }
            else
            {
                ionColor = MetaDrawSettings.ProductTypeToColor[matchedIon.NeutralTheoreticalProduct.ProductType];
            }

            int i = Scan.MassSpectrum.GetClosestPeakIndex(matchedIon.NeutralTheoreticalProduct.NeutralMass.ToMz(matchedIon.Charge));
            double mz = Scan.MassSpectrum.XArray[i];
            double intensity = Scan.MassSpectrum.YArray[i];

            if (useLiteralPassedValues)
            {
                mz = matchedIon.Mz;
                intensity = matchedIon.Intensity;
            }

            // peak annotation
            string prefix = "";
            if (SpectrumMatch.BetaPeptideBaseSequence != null)
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

            string peakAnnotationText = prefix + matchedIon.NeutralTheoreticalProduct.Annotation;

            if (matchedIon.NeutralTheoreticalProduct.NeutralLoss != 0)
            {
                peakAnnotationText += "-" + matchedIon.NeutralTheoreticalProduct.NeutralLoss.ToString("F2");
            }

            if (MetaDrawSettings.AnnotateCharges)
            {
                peakAnnotationText += "+" + matchedIon.Charge;
            }

            if (MetaDrawSettings.AnnotateMzValues)
            {
                peakAnnotationText += " (" + matchedIon.Mz.ToString("F3") + ")";
            }

            var peakAnnotation = new TextAnnotation();
            peakAnnotation.Font = "Arial";
            peakAnnotation.FontSize = MetaDrawSettings.AnnotatedFontSize;
            peakAnnotation.FontWeight = MetaDrawSettings.AnnotationBold ? OxyPlot.FontWeights.Bold : 2.0;
            peakAnnotation.TextColor = ionColor;
            peakAnnotation.StrokeThickness = 0;
            peakAnnotation.Text = peakAnnotationText;
            peakAnnotation.TextPosition = new DataPoint(mz, intensity);
            peakAnnotation.TextVerticalAlignment = intensity < 0 ? OxyPlot.VerticalAlignment.Top : OxyPlot.VerticalAlignment.Bottom;
            peakAnnotation.TextHorizontalAlignment = OxyPlot.HorizontalAlignment.Center;

            if (!MetaDrawSettings.DisplayIonAnnotations)
            {
                peakAnnotation.Text = string.Empty;
            }

            DrawPeak(mz, intensity, MetaDrawSettings.StrokeThicknessAnnotated, ionColor, peakAnnotation);
        }

        protected void AnnotateProperties()
        {
            StringBuilder text = new StringBuilder();
            if (MetaDrawSettings.SpectrumDescription["Precursor Charge: "])
            {
                text.Append("Precursor Charge: ");
                text.Append(SpectrumMatch.PrecursorCharge);
                text.Append("\r\n");
            }
            if (MetaDrawSettings.SpectrumDescription["Precursor Mass: "])
            {
                text.Append("Precursor Mass: ");
                text.Append(SpectrumMatch.PrecursorMass.ToString("F3"));
                text.Append("\r\n");
            }
            if (MetaDrawSettings.SpectrumDescription["Theoretical Mass: "])
            {
                text.Append("Theoretical Mass: ");
                text.Append(double.TryParse(SpectrumMatch.PeptideMonoMass, NumberStyles.Any, CultureInfo.InvariantCulture, out var monoMass) ? monoMass.ToString("F3") : SpectrumMatch.PeptideMonoMass);
                text.Append("\r\n");
            }
            if (MetaDrawSettings.SpectrumDescription["Protein Accession: "])
            {
                text.Append("Protein Accession: ");
                text.Append(SpectrumMatch.ProteinAccession);
                text.Append("\r\n");
            }
            if (SpectrumMatch.ProteinName != null && MetaDrawSettings.SpectrumDescription["Protein: "])
            {

                text.Append("Protein: ");
                text.Append(SpectrumMatch.ProteinName.Length > 20 ? SpectrumMatch.ProteinName.Substring(0, 18) + "..." : SpectrumMatch.ProteinName);
                text.Append("\r\n");

            }
            if (MetaDrawSettings.SpectrumDescription["Decoy/Contaminant/Target: "])
            {
                text.Append("Decoy/Contaminant/Target: ");
                text.Append(SpectrumMatch.DecoyContamTarget);
                text.Append("\r\n");
            }
            if (MetaDrawSettings.SpectrumDescription["Sequence Length: "])
            {
                text.Append("Sequence Length: ");
                text.Append(SpectrumMatch.BaseSeq.Length.ToString("F3").Split('.')[0]);
                text.Append("\r\n");
            }
            if (MetaDrawSettings.SpectrumDescription["ProForma Level: "])
            {
                text.Append("ProForma Level: ");
                text.Append(SpectrumMatch.AmbiguityLevel);
                text.Append("\r\n");
            }
            if (MetaDrawSettings.SpectrumDescription["Score: "])
            {
                text.Append("Score: ");
                text.Append(SpectrumMatch.Score.ToString("F3"));
                text.Append("\r\n");
            }
            if (MetaDrawSettings.SpectrumDescription["Q-Value: "])
            {
                text.Append("Q-Value: ");
                text.Append(SpectrumMatch.QValue.ToString("F3"));
                text.Append("\r\n");
            }
            if (MetaDrawSettings.SpectrumDescription["PEP: "])
            {
                text.Append("PEP: ");
                text.Append(SpectrumMatch.PEP);
                text.Append("\r\n");
            }
            if (MetaDrawSettings.SpectrumDescription["PEP Q-Value: "])
            {
                text.Append("PEP Q-Value: ");
                text.Append(SpectrumMatch.PEP_QValue);
                text.Append("\r\n");
            }

            var annotation = new PlotTextAnnotation()
            {
                Text = text.ToString(),
                XPosition = PlotTextAnnotation.RelativeToX.Right,
                YPosition = PlotTextAnnotation.RelativeToY.Top,
                FontWeight = 550,
                X = -155,
                Font = "Arial",
                FontSize = 10,
                TextColor = OxyColors.Black,
            };

            this.Model.Annotations.Add(annotation);
        }

        protected void DrawPeak(double mz, double intensity, double strokeWidth, OxyColor color, TextAnnotation annotation)
        {
            // peak line
            var line = new LineSeries();
            line.Color = color;
            line.StrokeThickness = strokeWidth;
            line.Points.Add(new DataPoint(mz, 0));
            line.Points.Add(new DataPoint(mz, intensity));

            if (annotation != null)
            {
                this.Model.Annotations.Add(annotation);
            }

            this.Model.Series.Add(line);
        }

        protected void ZoomAxes(IEnumerable<MatchedFragmentIon> matchedFragmentIons, double yZoom = 1.2)
        {
            double highestAnnotatedIntensity = 0;
            double highestAnnotatedMz = double.MinValue;
            double lowestAnnotatedMz = double.MaxValue;

            foreach (var ion in matchedFragmentIons)
            {
                double mz = ion.NeutralTheoreticalProduct.NeutralMass.ToMz(ion.Charge);
                int i = Scan.MassSpectrum.GetClosestPeakIndex(mz);
                double intensity = Scan.MassSpectrum.YArray[i];

                if (intensity > highestAnnotatedIntensity)
                {
                    highestAnnotatedIntensity = intensity;
                }

                if (highestAnnotatedMz < mz)
                {
                    highestAnnotatedMz = mz;
                }

                if (mz < lowestAnnotatedMz)
                {
                    lowestAnnotatedMz = mz;
                }
            }

            if (highestAnnotatedIntensity > 0)
            {
                this.Model.Axes[1].Zoom(0, highestAnnotatedIntensity * yZoom);
            }

            if (highestAnnotatedMz > double.MinValue && lowestAnnotatedMz < double.MaxValue)
            {
                this.Model.Axes[0].Zoom(lowestAnnotatedMz - 100, highestAnnotatedMz + 100);
            }
        }

    }

    //TODO: move this to mzLib (https://github.com/smith-chem-wisc/mzLib/blob/master/mzPlot/Annotations/PlotTextAnnotation.cs)
    public class PlotTextAnnotation : Annotation
    {
        public enum RelativeToX { Left, Center, Right }
        public enum RelativeToY { Top, Center, Bottom }

        /// <summary>
        /// Creates a "floating" text annotation, i.e., the annotation's x,y position is coupled to the chart area and not a data point.
        /// </summary>
        public PlotTextAnnotation()
        {

        }

        public string Text { get; set; }

        /// <summary>
        /// The x-position of the annotation, in pixels. Relative to the left of the chart (higher X is more to the right). Can be negative.
        /// </summary>
        public double X { get; set; }

        /// <summary>
        /// The y-position of the annotation, in pixels. Relative to the top of the chart (higher Y is lower on the chart). Can be negative.
        /// </summary>
        public double Y { get; set; }

        public RelativeToX XPosition { get; set; } = RelativeToX.Left;
        public RelativeToY YPosition { get; set; } = RelativeToY.Top;

        public override void Render(IRenderContext rc)
        {
            base.Render(rc);
            double pX = 0;
            double pY = 0;

            switch (XPosition)
            {
                case RelativeToX.Left: pX = PlotModel.PlotArea.Left; break;
                case RelativeToX.Center: pX = PlotModel.PlotArea.Center.X; break;
                case RelativeToX.Right: pX = PlotModel.PlotArea.Right; break;
            }

            switch (YPosition)
            {
                case RelativeToY.Top: pY = PlotModel.PlotArea.Top; break;
                case RelativeToY.Center: pY = PlotModel.PlotArea.Center.Y; break;
                case RelativeToY.Bottom: pY = PlotModel.PlotArea.Bottom; break;
            }

            pX += X;
            pY += Y;

            rc.DrawMultilineText(new ScreenPoint(pX, pY), Text, TextColor, Font, FontSize, FontWeight);
        }


    }
}