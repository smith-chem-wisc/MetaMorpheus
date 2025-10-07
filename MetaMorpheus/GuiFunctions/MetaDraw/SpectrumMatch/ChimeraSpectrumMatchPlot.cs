﻿using OxyPlot;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using OxyPlot.Wpf;
using Point = System.Windows.Point;
using Canvas = System.Windows.Controls.Canvas;
using System;
using GuiFunctions.MetaDraw;
using Easy.Common.Extensions;
using Omics.Fragmentation;

namespace GuiFunctions
{
    public class ChimeraSpectrumMatchPlot : SpectrumMatchPlot
    {
        public ChimeraSpectrumMatchPlot(PlotView plotView, ChimeraGroupViewModel chimeraGroupVm, double mzMax = double.MaxValue) : base(plotView, null, chimeraGroupVm.Ms2Scan)
        {
            var matchedIonsByColor = chimeraGroupVm.MatchedFragmentIonsByColor;
            int totalCount = 0;
            foreach (var group in matchedIonsByColor.Values)
                totalCount += group.Count;
            var matchedIons = new List<MatchedFragmentIon>(totalCount);
            foreach (var group in matchedIonsByColor.Values)
                for (int i = 0; i < group.Count; i++)
                    matchedIons.Add(group[i].Item1);
            MatchedFragmentIons = matchedIons;

            ZoomAxes();
            AnnotateMatchedIons(chimeraGroupVm);
            if (Math.Abs(mzMax - double.MaxValue) > 0.001)
            {
                Model.Axes[0].Maximum = mzMax;
            }

            RefreshChart();
        }

        /// <summary>
        /// Annotates the matched ions based upon the protein of origin, and the unique proteoform ID's
        /// </summary>
        private void AnnotateMatchedIons(ChimeraGroupViewModel chimeraGroupVm)
        {
            foreach (var ionGroup in chimeraGroupVm.MatchedFragmentIonsByColor)
            {
                var color = ionGroup.Key;
                var ions = ionGroup.Value;
                for (int i = 0; i < ions.Count; i++)
                {
                    var tuple = ions[i];
                    AnnotatePeak(tuple.Item1, false, false, color);
                }
            }
        }

        public void ExportPlot(string path, Canvas legend = null, Point? legendPoint = null,  double width = 700, double height = 370)
        {
            width = width > 0 ? width : 700;
            height = height > 0 ? height : 300;
            var tempModelPath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(path), "temp." + MetaDrawSettings.ExportType);
            var tempLegendPngPath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(path), "legend.png");
            List<System.Drawing.Bitmap> bitmaps = new();
            List<Point> points = new();
            var dpiScale = MetaDrawSettings.CanvasPdfExportDpi / 96.0;

            // export model as png and load as bitmap
            ExportToPng(tempModelPath, (int)width, (int)height);
            bitmaps.Add(new System.Drawing.Bitmap(tempModelPath));
            points.Add(new Point(0, 0));

            // Render legend as bitmap and export as png if used
            System.Drawing.Bitmap ptmLegendBitmap = null;

            if (legend != null && MetaDrawSettings.ShowLegend && legendPoint.HasValue)
            {
                RenderTargetBitmap legendRenderBitmap = new((int)(dpiScale * legend.ActualWidth), (int)(dpiScale * legend.ActualHeight),
                MetaDrawSettings.CanvasPdfExportDpi, MetaDrawSettings.CanvasPdfExportDpi, PixelFormats.Pbgra32);
                legendRenderBitmap.Render(legend);
                var legendEncoder = new PngBitmapEncoder();
                legendEncoder.Frames.Add(BitmapFrame.Create(legendRenderBitmap));
                using (var file = File.Create(tempLegendPngPath))
                {
                    legendEncoder.Save(file);
                }

                System.Drawing.Bitmap tempLegendBitmap = new(tempLegendPngPath);
                ptmLegendBitmap = new System.Drawing.Bitmap(tempLegendBitmap, new System.Drawing.Size((int)legend.ActualWidth, (int)legend.ActualHeight));
                bitmaps.Add(ptmLegendBitmap);
                points.Add(legendPoint.Value);
                tempLegendBitmap.Dispose();
            }

            // combine the bitmaps
            var combinedBitmaps = MetaDrawLogic.CombineBitmap(bitmaps, points, true);
            bitmaps.ForEach(p => p.Dispose());
            File.Delete(tempModelPath);
            File.Delete(tempLegendPngPath);
            ExportPlot(path, combinedBitmaps, width, height);
        }
    }
}
