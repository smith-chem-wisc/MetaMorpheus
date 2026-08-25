using EngineLayer;
using MassSpectrometry;
using System.Collections.Generic;
using System.IO;
using System.Reflection;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using Readers;
using System;

namespace GuiFunctions
{
    /// <summary>
    /// Class for the peptide spectrum match plot within the metadraw window
    /// </summary>
    public class PeptideSpectrumMatchPlot : SpectrumMatchPlot
    {
        public PeptideSpectrumMatchPlot(OxyPlot.Wpf.PlotView plotView, SpectrumMatchFromTsv sm, MsDataScan scan,
            List<MatchedFragmentIon> matchedFragmentIons, bool annotateProperties = true, LibrarySpectrum librarySpectrum = null, bool stationarySequence = false)
            : base(plotView, sm, scan, matchedFragmentIons)
        {
            if (annotateProperties)
            {
                if (librarySpectrum != null)
                {
                    AnnotateProperties(librarySpectrum);
                }
                else
                {
                    AnnotateProperties();
                }
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
            // Use unique temp file names to ensure thread safety
            string tempModelPath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(path), Path.GetRandomFileName() + "." + MetaDrawSettings.ExportType);
            string tempStationarySequencePngPath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(path), Path.GetRandomFileName() + "_annotation.png");
            string tempPtmLegendPngPath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(path), Path.GetRandomFileName() + "_legend.png");
            List<System.Drawing.Bitmap> bitmaps = new();
            List<Point> points = new();

            // scales for desired DPI
            double dpiScale = MetaDrawSettings.CanvasPdfExportDpi / 96.0;

            // render stationary sequence as bitmap and export as png
            int stationarySequenceWidth = MetaDrawLogic.GetCanvasDimension(stationarySequence.Width, stationarySequence.ActualWidth) + 30;
            int stationarySequenceHeight = MetaDrawLogic.GetCanvasDimension(stationarySequence.Height, stationarySequence.ActualHeight) + 30;
            Size stationarySequenceSize = new Size(stationarySequenceWidth, stationarySequenceHeight);
            double originalStationarySequenceWidth = stationarySequence.Width;
            double originalStationarySequenceHeight = stationarySequence.Height;
            RenderTargetBitmap renderStationaryBitmap;
            Vector stationarySequenceLocationVector;
            try
            {
                stationarySequence.Width = stationarySequenceWidth;
                stationarySequence.Height = stationarySequenceHeight;
                stationarySequence.Measure(stationarySequenceSize);
                stationarySequence.Arrange(new Rect(stationarySequenceSize));

                renderStationaryBitmap = new RenderTargetBitmap((int)(dpiScale * stationarySequenceWidth), (int)(dpiScale * stationarySequenceHeight),
                                                      MetaDrawSettings.CanvasPdfExportDpi, MetaDrawSettings.CanvasPdfExportDpi, PixelFormats.Pbgra32);
                renderStationaryBitmap.Render(stationarySequence);
                stationarySequenceLocationVector = (Vector)stationarySequence.GetType().GetProperty("VisualOffset", BindingFlags.NonPublic | BindingFlags.Instance).GetValue(stationarySequence);
            }
            finally
            {
                stationarySequence.Width = originalStationarySequenceWidth;
                stationarySequence.Height = originalStationarySequenceHeight;
            }

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
            System.Drawing.Bitmap stationarySequenceBitmap = new System.Drawing.Bitmap(tempStatSequenceBitmap, new System.Drawing.Size(stationarySequenceWidth, stationarySequenceHeight));
            bitmaps.Add(stationarySequenceBitmap);

            // Magic Number +50 is used to offset the stationary sequence to the right of the model plot. 
            Point stationarySequencePoint = new Point(stationarySequenceLocationVector.X + 50, stationarySequenceLocationVector.Y);
            points.Add(stationarySequencePoint);

            // render ptm legend as bitmap and export as png if used
            System.Drawing.Bitmap ptmLegendBitmap = null;
            Point ptmLegendPoint;
            if (ptmLegend != null && MetaDrawSettings.ShowLegend)
            {
                // Saving Canvas as a usable Png
                int ptmLegendWidth = MetaDrawLogic.GetCanvasDimension(ptmLegend.ActualWidth, 0);
                int ptmLegendHeight = MetaDrawLogic.GetCanvasDimension(ptmLegend.ActualHeight, 0);
                RenderTargetBitmap ptmLegendRenderBitmap = new((int)(dpiScale * ptmLegendWidth), (int)(dpiScale * ptmLegendHeight),
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
                ptmLegendBitmap = new System.Drawing.Bitmap(tempPtmLegendBitmap, new System.Drawing.Size(ptmLegendWidth, ptmLegendHeight));
                bitmaps.Add(ptmLegendBitmap);

                // Magic Number -14 and -20 is used to offset the ptm legend to the left and above the bottom left corner of the model plot.
                ptmLegendPoint = new Point(ptmLegendLocationVector.X - 14, ptmLegendLocationVector.Y - 20);
                points.Add(ptmLegendPoint);
                tempPtmLegendBitmap.Dispose();
            }

            // combine the bitmaps
            System.Drawing.Bitmap combinedBitmaps = MetaDrawLogic.CombineBitmap(bitmaps, points);
            tempStatSequenceBitmap.Dispose();
            File.Delete(tempModelPath);
            File.Delete(tempStationarySequencePngPath);
            File.Delete(tempPtmLegendPngPath);
            ExportPlot(path, combinedBitmaps, width, height);
        }
    }
}
