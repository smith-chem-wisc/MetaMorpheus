using Chemistry;
using iText.IO.Image;
using iText.Kernel.Pdf;
using MassSpectrometry;
using mzPlot;
using OxyPlot;
using OxyPlot.Annotations;
using OxyPlot.Axes;
using OxyPlot.Series;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;

namespace EngineLayer
{
    public class PeptideSpectrumMatchPlot : Plot
    {
        public MsDataScan Scan { get; protected set; }
        protected PsmFromTsv SpectrumMatch;
        protected Canvas SequenceDrawingCanvas;

        public PeptideSpectrumMatchPlot(OxyPlot.Wpf.PlotView plotView, Canvas sequenceDrawingCanvas, PsmFromTsv psm, MsDataScan scan,
            List<MatchedFragmentIon> matchedFragmentIons, bool annotateProperties = true) : base(plotView)
        {
            Model.Title = string.Empty;
            Model.Subtitle = string.Empty;
            SpectrumMatch = psm;
            Scan = scan;
            SequenceDrawingCanvas = sequenceDrawingCanvas;
            SequenceDrawingCanvas.Height = 60;
            sequenceDrawingCanvas.Width = 600;
            
            ClearCanvas(SequenceDrawingCanvas);
            DrawSpectrum();
            AnnotateBaseSequence(psm.BaseSeq, psm.FullSequence, 10, matchedFragmentIons);
            AnnotateMatchedIons(isBetaPeptide: false, matchedFragmentIons);

            if (annotateProperties)
            {
                AnnotateProperties();
            }

            ZoomAxes(matchedFragmentIons);
            RefreshChart();
        }

        public new void ExportToPdf(string path, double width = 700, double height = 370)
        {
            // exports spectrum annotation w/o base seq annotation
            string tempPdfPath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(path), "temp.pdf");
            string tempPngPath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(path), "annotation.png");
            base.ExportToPdf(tempPdfPath, width, height);

            // scales for desired DPI
            double dpiScale = MetaDrawSettings.CanvasPdfExportDpi / 96.0;

            // save base seq as PNG
            SequenceDrawingCanvas.Measure(new Size((int)SequenceDrawingCanvas.Width, (int)SequenceDrawingCanvas.Height));
            SequenceDrawingCanvas.Arrange(new Rect(new Size((int)SequenceDrawingCanvas.Width, (int)SequenceDrawingCanvas.Height)));

            RenderTargetBitmap renderBitmap = new RenderTargetBitmap((int)(dpiScale * SequenceDrawingCanvas.Width), (int)(dpiScale * SequenceDrawingCanvas.Height), 
                MetaDrawSettings.CanvasPdfExportDpi, MetaDrawSettings.CanvasPdfExportDpi, PixelFormats.Pbgra32);

            renderBitmap.Render(SequenceDrawingCanvas);
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
            img.SetMarginLeft((float)(-1.0 * SequenceDrawingCanvas.Margin.Left) + 10);
            img.SetMarginTop(-30);
            img.ScaleToFit((float)SequenceDrawingCanvas.Width, (float)SequenceDrawingCanvas.Height);

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
                AxisTitleDistance = 10
            });

            // draw all peaks in the scan
            for (int i = 0; i < Scan.MassSpectrum.XArray.Length; i++)
            {
                double mz = Scan.MassSpectrum.XArray[i];
                double intensity = Scan.MassSpectrum.YArray[i];

                DrawPeak(mz, intensity, MetaDrawSettings.StrokeThicknessUnannotated, MetaDrawSettings.UnannotatedPeakColor, null);
            }
        }

        protected void AnnotateMatchedIons(bool isBetaPeptide, List<MatchedFragmentIon> matchedFragmentIons)
        {
            foreach (MatchedFragmentIon matchedIon in matchedFragmentIons)
            {
                AnnotatePeak(matchedIon, isBetaPeptide);
            }
        }

        protected void AnnotateBaseSequence(string baseSequence, string fullSequence, int yLoc, List<MatchedFragmentIon> matchedFragmentIons)
        {
            // don't draw ambiguous sequences
            if (SpectrumMatch.FullSequence.Contains("|"))
            {
                return;
            }

            // draw base sequence
            double canvasWidth = SequenceDrawingCanvas.Width;
            for (int r = 0; r < baseSequence.Length; r++)
            {
                double x = r * MetaDrawSettings.AnnotatedSequenceTextSpacing + 10;
                DrawText(SequenceDrawingCanvas, new Point(x, yLoc), baseSequence[r].ToString(), Brushes.Black);

                canvasWidth = x + 30;
            }
            SequenceDrawingCanvas.Width = Math.Max(SequenceDrawingCanvas.Width, canvasWidth);

            // draw the fragment ion annotations on the base sequence
            foreach (var ion in matchedFragmentIons)
            {
                int residue = ion.NeutralTheoreticalProduct.AminoAcidPosition;
                string annotation = ion.NeutralTheoreticalProduct.ProductType + "" + ion.NeutralTheoreticalProduct.FragmentNumber;
                OxyColor oxycolor = SpectrumMatch.VariantCrossingIons.Contains(ion) ?
                    MetaDrawSettings.variantCrossColor : MetaDrawSettings.productTypeToColor[ion.NeutralTheoreticalProduct.ProductType];
                Color color = Color.FromArgb(oxycolor.A, oxycolor.R, oxycolor.G, oxycolor.B);

                if (ion.NeutralTheoreticalProduct.NeutralLoss != 0)
                {
                    annotation += "-" + ion.NeutralTheoreticalProduct.NeutralLoss;
                }

                double x = residue * MetaDrawSettings.AnnotatedSequenceTextSpacing + 11;
                double y = yLoc + MetaDrawSettings.productTypeToYOffset[ion.NeutralTheoreticalProduct.ProductType];

                if (ion.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.C)
                {
                    DrawCTermIon(SequenceDrawingCanvas, new Point(x, y), color, annotation);
                }
                else if (ion.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.N)
                {
                    DrawNTermIon(SequenceDrawingCanvas, new Point(x, y), color, annotation);
                }
                // don't draw diagnostic ions, precursor ions, etc
            }

            AnnotateModifications(fullSequence, yLoc);
        }

        protected void AnnotateModifications(string fullSequence, int yLoc)
        {
            var peptide = new PeptideWithSetModifications(fullSequence, GlobalVariables.AllModsKnownDictionary);

            // read glycans if applicable
            List<Tuple<int, string, double>> localGlycans = null;
            if (SpectrumMatch.GlycanLocalizationLevel != null)
            {
                localGlycans = PsmFromTsv.ReadLocalizedGlycan(SpectrumMatch.LocalizedGlycan);
            }

            // annotate mods
            foreach (var mod in peptide.AllModsOneIsNterminus)
            {
                double xLocation = (mod.Key - 1) * MetaDrawSettings.AnnotatedSequenceTextSpacing - 12;
                double yLocation = yLoc + 2;

                if (mod.Value.ModificationType == "O-Glycosylation")
                {
                    if (localGlycans.Where(p => p.Item1 + 1 == mod.Key).Count() > 0)
                    {
                        DrawCircle(SequenceDrawingCanvas, new Point(xLocation, yLocation), MetaDrawSettings.ModificationAnnotationColor);
                    }
                    else
                    {
                        DrawCircle(SequenceDrawingCanvas, new Point(xLocation, yLocation), Brushes.Gray);
                    }
                }
                else
                {
                    DrawCircle(SequenceDrawingCanvas, new Point(xLocation, yLocation), MetaDrawSettings.ModificationAnnotationColor);
                }
            }
        }

        protected void AnnotatePeak(MatchedFragmentIon matchedIon, bool isBetaPeptide)
        {
            OxyColor ionColor;

            if (SpectrumMatch.VariantCrossingIons.Contains(matchedIon))
            {
                ionColor = MetaDrawSettings.variantCrossColor;
            }
            else
            {
                if (isBetaPeptide)
                {
                    ionColor = MetaDrawSettings.betaProductTypeToColor[matchedIon.NeutralTheoreticalProduct.ProductType];
                }
                else
                {
                    ionColor = MetaDrawSettings.productTypeToColor[matchedIon.NeutralTheoreticalProduct.ProductType];
                }
            }

            int i = Scan.MassSpectrum.GetClosestPeakIndex(matchedIon.NeutralTheoreticalProduct.NeutralMass.ToMz(matchedIon.Charge));
            double mz = matchedIon.Mz;
            double intensity = Scan.MassSpectrum.YArray[i];

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

            string productType = matchedIon.NeutralTheoreticalProduct.ProductType.ToString()
                //.Replace("star", "*", StringComparison.OrdinalIgnoreCase)
                //.Replace("degree", "°", StringComparison.OrdinalIgnoreCase)
                //.Replace("dot", "·", StringComparison.OrdinalIgnoreCase)
                ;
            string productNumber = matchedIon.NeutralTheoreticalProduct.FragmentNumber.ToString();
            string peakAnnotationText = prefix + productType + productNumber;

            if (matchedIon.NeutralTheoreticalProduct.NeutralLoss != 0)
            {
                peakAnnotationText += "-" + matchedIon.NeutralTheoreticalProduct.NeutralLoss.ToString("F2");
            }

            if (MetaDrawSettings.ShowAnnotationCharges)
            {
                peakAnnotationText += "+" + matchedIon.Charge;
            }

            if (MetaDrawSettings.ShowMzValues)
            {
                peakAnnotationText += " (" + matchedIon.Mz.ToString("F3") + ")";
            }

            var peakAnnotation = new TextAnnotation();
            peakAnnotation.Font = "Arial";
            peakAnnotation.FontSize = MetaDrawSettings.AnnotatedFontSize;
            peakAnnotation.FontWeight = MetaDrawSettings.BoldText ? OxyPlot.FontWeights.Bold : 2.0;
            peakAnnotation.TextColor = ionColor;
            peakAnnotation.StrokeThickness = 0;
            peakAnnotation.Text = peakAnnotationText;
            peakAnnotation.TextPosition = new DataPoint(mz, intensity);
            peakAnnotation.TextHorizontalAlignment = OxyPlot.HorizontalAlignment.Center;

            DrawPeak(mz, intensity, MetaDrawSettings.StrokeThicknessAnnotated, ionColor, peakAnnotation);
        }

        protected void AnnotateProperties()
        {
            StringBuilder text = new StringBuilder();
            text.Append("Precursor Charge: ");
            text.Append(SpectrumMatch.PrecursorCharge);
            text.Append("\r\n");

            text.Append("Precursor Mass: ");
            text.Append(SpectrumMatch.PrecursorMass.ToString("F3"));
            text.Append("\r\n");

            text.Append("Theoretical Mass: ");
            text.Append(double.TryParse(SpectrumMatch.PeptideMonoMass, out var monoMass) ? monoMass.ToString("F3") : SpectrumMatch.PeptideMonoMass);
            text.Append("\r\n");

            text.Append("Score: ");
            text.Append(SpectrumMatch.Score.ToString("F3"));
            text.Append("\r\n");

            text.Append("Protein Accession: ");
            text.Append(SpectrumMatch.ProteinAccession);
            text.Append("\r\n");

            if (SpectrumMatch.ProteinName != null)
            {
                text.Append("Protein: ");
                text.Append(SpectrumMatch.ProteinName.Length > 20 ? SpectrumMatch.ProteinName.Substring(0, 18) + "..." : SpectrumMatch.ProteinName);
                text.Append("\r\n");
            }

            text.Append("Decoy/Contaminant/Target: ");
            text.Append(SpectrumMatch.DecoyContamTarget);
            text.Append("\r\n");

            text.Append("Q-Value: ");
            text.Append(SpectrumMatch.QValue.ToString("F3"));
            text.Append("\r\n");

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

        protected void ZoomAxes(IEnumerable<MatchedFragmentIon> matchedFragmentIons, double yZoom  = 1.2)
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

        /// <summary>
        /// Draw the line seperator @ top
        /// </summary>
        private static void DrawCTermIon(Canvas cav, Point topLoc, Color clr, string footnote)
        {
            double x = topLoc.X, y = topLoc.Y;
            Polyline bot = new Polyline();
            bot.Points = new PointCollection() { new Point(x + 10, y + 10), new Point(x, y + 10), new Point(x, y + 24) };
            bot.Stroke = new SolidColorBrush(clr);
            bot.StrokeThickness = 1;
            cav.Children.Add(bot);
            Canvas.SetZIndex(bot, 1); //on top of any other things in canvas
        }

        /// <summary>
        /// Draw the line seperator @ bottom
        /// </summary>
        private static void DrawNTermIon(Canvas cav, Point botLoc, Color clr, string footnote)
        {
            double x = botLoc.X, y = botLoc.Y;
            Polyline bot = new Polyline();
            bot.Points = new PointCollection() { new Point(x - 10, y - 10), new Point(x, y - 10), new Point(x, y - 24) };
            bot.Stroke = new SolidColorBrush(clr);
            bot.StrokeThickness = 1;
            Canvas.SetZIndex(bot, 1); //on top of any other things in canvas
            cav.Children.Add(bot);
        }

        /// <summary>
        /// Create text blocks on canvas
        /// </summary>
        private static void DrawText(Canvas cav, Point loc, string txt, Brush clr)
        {
            TextBlock tb = new TextBlock();
            tb.Foreground = clr;
            tb.Text = txt;
            tb.Height = 30;
            tb.FontSize = 25;
            tb.FontWeight = System.Windows.FontWeights.Bold;
            tb.FontFamily = new FontFamily("Arial");
            tb.TextAlignment = TextAlignment.Center;
            tb.HorizontalAlignment = System.Windows.HorizontalAlignment.Center;
            tb.Width = 24; // W (tryptophan) seems to be widest letter, make sure it fits if you're editing this

            Canvas.SetTop(tb, loc.Y);
            Canvas.SetLeft(tb, loc.X);
            Panel.SetZIndex(tb, 2); //lower priority
            cav.Children.Add(tb);
            cav.UpdateLayout();
        }

        /// <summary>
        /// Draws a circle
        /// </summary>
        private static void DrawCircle(Canvas cav, Point loc, SolidColorBrush clr)
        {
            Ellipse circle = new Ellipse()
            {
                Width = 24,
                Height = 24,
                Stroke = clr,
                StrokeThickness = 1,
                Fill = clr,
                Opacity = 0.7
            };
            Canvas.SetLeft(circle, loc.X);
            Canvas.SetTop(circle, loc.Y);
            Panel.SetZIndex(circle, 1);
            cav.Children.Add(circle);
        }

        /// <summary>
        /// Clear canvas board
        /// </summary>
        private static void ClearCanvas(Canvas cav)
        {
            cav.Children.Clear();
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