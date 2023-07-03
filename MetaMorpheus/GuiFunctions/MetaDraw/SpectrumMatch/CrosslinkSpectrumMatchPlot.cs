using System.Collections.Generic;
using System.IO;
using EngineLayer;
using MassSpectrometry;
using System.Linq;
using System.Reflection;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Media.Imaging;

namespace GuiFunctions
{
    /// <summary>
    /// Class for displaying CrossLinked spectra in MetaDraw
    /// </summary>
    public class CrosslinkSpectrumMatchPlot : SpectrumMatchPlot
    {
        public CrosslinkSpectrumMatchPlot(OxyPlot.Wpf.PlotView plotView, PsmFromTsv csm, MsDataScan scan, Canvas stationaryCanvas,
            bool annotateProperties = true, LibrarySpectrum librarySpectrum = null)
            : base(plotView, csm, scan)
        {

            if (annotateProperties)
            {
                AnnotateProperties(librarySpectrum);
            }

            // annotate beta peptide matched ions
            AnnotateMatchedIons(isBetaPeptide: true, csm.BetaPeptideMatchedIons);

            ZoomAxes(csm.MatchedIons.Concat(csm.BetaPeptideMatchedIons), yZoom: 1.5);
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
            ExportToPng(tempModelPath, (int)width, (int)height);
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
            File.Delete(tempModelPath);
            File.Delete(tempStationarySequencePngPath);
            File.Delete(tempPtmLegendPngPath);
            base.ExportPlot(path, combinedBitmaps, width, height);
        }
    }
}