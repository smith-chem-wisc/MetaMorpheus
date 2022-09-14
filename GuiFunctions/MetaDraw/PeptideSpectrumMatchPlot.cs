using Chemistry;
using EngineLayer;
using iText.IO.Image;
using iText.Kernel.Pdf;
using iText.Layout;
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
using System.Reflection;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using iText.Kernel.Geom;
using Point = System.Windows.Point;
using Vector = System.Windows.Vector;
using Canvas = System.Windows.Controls.Canvas;

namespace GuiFunctions
{
    /// <summary>
    /// Class for the peptide spectrum match plot within the metadraw window
    /// </summary>
    public class PeptideSpectrumMatchPlot : SpectrumMatchPlot
    {
        public PeptideSpectrumMatchPlot(OxyPlot.Wpf.PlotView plotView, PsmFromTsv psm, MsDataScan scan,
            List<MatchedFragmentIon> matchedFragmentIons, bool annotateProperties = true, LibrarySpectrum librarySpectrum = null, bool stationarySequence = false) 
            : base(plotView, psm, scan)
        {
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

        public void ExportPlot(string path, Canvas stationarySequence, Canvas ptmLegend = null, Vector ptmLegendLocationVector = new(), double width = 700, double height = 370)
        {
            width = width > 0 ? width : 700;
            height = height > 0 ? height : 300;
            string tempModelPath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(path), "temp." + MetaDrawSettings.ExportType);
            string tempStationarySequencePngPath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(path), "annotation.png");
            string tempPtmLegendPngPath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(path), "legend.png");
            List<System.Drawing.Bitmap> bitmaps = new();
            List<Point> points = new();

            // scales for desired DPI
            double dpiScale = MetaDrawSettings.CanvasPdfExportDpi / 96.0;

            // render stationary sequence as bitmap and export as png
            stationarySequence.Height += 30;
            stationarySequence.Width += 30;
            Size stationarySequenceSize = new Size((int)stationarySequence.Width, (int)stationarySequence.Height);
            stationarySequence.Measure(stationarySequenceSize);
            stationarySequence.Arrange(new Rect(stationarySequenceSize));

            RenderTargetBitmap renderStationaryBitmap = new RenderTargetBitmap((int)(dpiScale * stationarySequence.Width), (int)(dpiScale * stationarySequence.Height),
                                                  MetaDrawSettings.CanvasPdfExportDpi, MetaDrawSettings.CanvasPdfExportDpi, PixelFormats.Pbgra32);
            renderStationaryBitmap.Render(stationarySequence);

            PngBitmapEncoder encoder = new PngBitmapEncoder();
            encoder.Frames.Add(BitmapFrame.Create(renderStationaryBitmap));
            using (FileStream file = File.Create(tempStationarySequencePngPath))
            {
                encoder.Save(file);
            }

            // export model as png and load both stationary and model as bitmap
            base.ExportToPng(tempModelPath, (int)width, (int)height);
            bitmaps.Add(new System.Drawing.Bitmap(tempModelPath));
            points.Add(new Point(0, 0));

            var tempStatSequenceBitmap = new System.Drawing.Bitmap(tempStationarySequencePngPath);
            System.Drawing.Bitmap stationarySequenceBitmap = new System.Drawing.Bitmap(tempStatSequenceBitmap, new System.Drawing.Size((int)stationarySequence.Width, (int)stationarySequence.Height));
            bitmaps.Add(stationarySequenceBitmap);
            var stationarySequenceLocationVector = (Vector)stationarySequence.GetType().GetProperty("VisualOffset", BindingFlags.NonPublic | BindingFlags.Instance).GetValue(stationarySequence);
            Point stationarySequencePoint = new Point(stationarySequenceLocationVector.X, stationarySequenceLocationVector.Y);
            points.Add(stationarySequencePoint);

            // render ptm legend as bitmap and export as png if used
            System.Drawing.Bitmap ptmLegendBitmap = null;
            Point ptmLegendPoint;
            if (ptmLegend != null && MetaDrawSettings.ShowLegend)
            {
                // Saving Canvas as a usable Png
                RenderTargetBitmap ptmLegendRenderBitmap = new((int)(dpiScale * ptmLegend.ActualWidth), (int)(dpiScale * ptmLegend.ActualHeight),
                         MetaDrawSettings.CanvasPdfExportDpi, MetaDrawSettings.CanvasPdfExportDpi, PixelFormats.Pbgra32);
                ptmLegendRenderBitmap.Render(ptmLegend);
                PngBitmapEncoder legendEncoder = new PngBitmapEncoder();
                legendEncoder.Frames.Add(BitmapFrame.Create(ptmLegendRenderBitmap));
                using (FileStream file = File.Create(tempPtmLegendPngPath))
                {
                    legendEncoder.Save(file);
                }

                // converting png to the final bitmap format
                System.Drawing.Bitmap tempPtmLegendBitmap = new(tempPtmLegendPngPath);
                ptmLegendBitmap = new System.Drawing.Bitmap(tempPtmLegendBitmap, new System.Drawing.Size((int)ptmLegend.ActualWidth, (int)ptmLegend.ActualHeight));
                bitmaps.Add(ptmLegendBitmap);
                ptmLegendPoint = new Point(ptmLegendLocationVector.X, ptmLegendLocationVector.Y);
                points.Add(ptmLegendPoint);
                tempPtmLegendBitmap.Dispose();
            }

            // combine the bitmaps
            System.Drawing.Bitmap combinedBitmaps = MetaDrawLogic.CombineBitmap(bitmaps, points);
            tempStatSequenceBitmap.Dispose();

            switch (MetaDrawSettings.ExportType)
            {
                case "Pdf":
                    string tempCombinedPath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(path), "tempCombined.png");
                    combinedBitmaps.Save(tempCombinedPath, System.Drawing.Imaging.ImageFormat.Png);

                    PdfDocument pdfDoc = new(new PdfWriter(path));
                    iText.Layout.Document document = new(pdfDoc, new iText.Kernel.Geom.PageSize((float)width - 30, (float)height - 30));

                    ImageData sequenceAndLegendImageData = ImageDataFactory.Create(tempCombinedPath);
                    iText.Layout.Element.Image sequenceAndPtmLegendImage = new(sequenceAndLegendImageData);
                    sequenceAndPtmLegendImage.SetMarginLeft(-30);
                    sequenceAndPtmLegendImage.SetMarginTop(-30);
                    sequenceAndPtmLegendImage.ScaleToFit((float)width, (float)height);
                    document.Add(sequenceAndPtmLegendImage);

                    pdfDoc.Close();
                    document.Close();
                    File.Delete(tempCombinedPath);
                    break;

                case "Png":
                    combinedBitmaps.Save(path, System.Drawing.Imaging.ImageFormat.Png);
                    break;

                case "Jpeg":
                    combinedBitmaps.Save(path, System.Drawing.Imaging.ImageFormat.Jpeg);
                    break;

                case "Tiff":
                    combinedBitmaps.Save(path, System.Drawing.Imaging.ImageFormat.Tiff);
                    break;

                case "Wmf":
                    combinedBitmaps.Save(path, System.Drawing.Imaging.ImageFormat.Wmf);
                    break;

                case "Bmp":
                    combinedBitmaps.Save(path, System.Drawing.Imaging.ImageFormat.Bmp);
                    break;
            }

            File.Delete(tempModelPath);
            File.Delete(tempStationarySequencePngPath);
            File.Delete(tempPtmLegendPngPath);
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
                if (SpectrumMatch.ProteinAccession.Length > 10)
                {
                    text.Append("\r\n   " + SpectrumMatch.ProteinAccession);
                }
                else
                    text.Append(SpectrumMatch.ProteinAccession);
                text.Append("\r\n");
            }
            if (SpectrumMatch.ProteinName != null && MetaDrawSettings.SpectrumDescription["Protein: "])
            {
                text.Append("Protein: ");
                if (SpectrumMatch.ProteinName.Length > 20)
                {
                    text.Append(SpectrumMatch.ProteinName.Substring(0, 20));
                    int length = SpectrumMatch.ProteinName.Length;
                    int remaining = length - 20;
                    for (int i = 20; i < SpectrumMatch.ProteinName.Length; i += 26)
                    {
                        if (remaining <= 26)
                            text.Append("\r\n   " + SpectrumMatch.ProteinName.Substring(i, remaining - 1));
                        else
                        {
                            text.Append("\r\n   " + SpectrumMatch.ProteinName.Substring(i, 26));
                            remaining -= 26;
                        }
                    }
                }
                else
                    text.Append(SpectrumMatch.ProteinName);
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
            if (MetaDrawSettings.SpectrumDescription["Spectral Angle: "])
            {
                text.Append("Spectral Angle: ");
                if (SpectrumMatch.SpectralAngle != null)
                    text.Append(SpectrumMatch.SpectralAngle.ToString());
                else
                    text.Append("N/A");
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
                text.Append(SpectrumMatch.PEP.ToString("F3"));
                text.Append("\r\n");
            }
            if (MetaDrawSettings.SpectrumDescription["PEP Q-Value: "])
            {
                text.Append("PEP Q-Value: ");
                text.Append(SpectrumMatch.PEP_QValue.ToString("F3"));
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