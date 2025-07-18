using OxyPlot;
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

namespace GuiFunctions
{
    public class ChimeraSpectrumMatchPlot : SpectrumMatchPlot
    {
        private static Queue<OxyColor> overflowColors;
        public static OxyColor MultipleProteinSharedColor;
        public static Dictionary<int, List<OxyColor>> ColorByProteinDictionary;
        public static Queue<OxyColor> OverflowColors => new(overflowColors);

        public ChimeraSpectrumMatchPlot(PlotView plotView, ChimeraGroupViewModel chimeraGroupVm, double mzMax = double.MaxValue) : base(plotView, null, chimeraGroupVm.Ms2Scan)
        {
            MatchedFragmentIons = chimeraGroupVm.MatchedFragmentIonsByColor.SelectMany(p => p.Value.Select(q => q.Item1)).ToList();
            AnnotateMatchedIons(chimeraGroupVm);
            if (Math.Abs(mzMax - double.MaxValue) > 0.001)
            {
                Model.Axes[0].Maximum = mzMax;
            }

            ZoomAxes();
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
                ionGroup.Value.ForEach(p => AnnotatePeak(p.Item1, false, false, color));
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

        /// <summary>
        /// Initializes the colors to be used by the Chimera Plotter
        /// </summary>
        static ChimeraSpectrumMatchPlot()
        {
            MultipleProteinSharedColor = OxyColors.Black;
            ColorByProteinDictionary = new()
            {
                {
                    0, [
                        OxyColors.Blue, OxyColors.MediumBlue, OxyColors.CornflowerBlue,
                        OxyColors.DarkBlue, OxyColors.CadetBlue, OxyColors.SteelBlue, OxyColors.DodgerBlue,
                        OxyColors.AliceBlue, OxyColors.DarkSlateBlue, OxyColors.DeepSkyBlue, OxyColors.DodgerBlue
                    ]
                },
                {
                    1, [
                        OxyColors.Red, OxyColors.IndianRed, OxyColors.PaleVioletRed,
                        OxyColors.LightCoral, OxyColors.Firebrick, OxyColors.Maroon, OxyColors.Tomato
                    ]
                },
                {
                    2, [
                        OxyColors.Green, OxyColors.MediumSpringGreen, OxyColors.LightGreen,
                        OxyColors.Linen, OxyColors.SpringGreen, OxyColors.Chartreuse, OxyColors.DarkSeaGreen
                    ]
                },
                {
                    3, [
                        OxyColors.Purple, OxyColors.MediumPurple, OxyColors.Violet,
                        OxyColors.Plum, OxyColors.Orchid, OxyColors.BlueViolet, OxyColors.Magenta
                    ]
                },
                {
                    4, [
                        OxyColors.Brown, OxyColors.SaddleBrown, OxyColors.Sienna, OxyColors.Chocolate,
                        OxyColors.SandyBrown, OxyColors.Chocolate, OxyColors.Peru, OxyColors.Tan
                    ]
                },
                {
                    5, [
                        OxyColors.Gold, OxyColors.DarkGoldenrod, OxyColors.Wheat, OxyColors.Goldenrod,
                        OxyColors.DarkKhaki, OxyColors.Khaki, OxyColors.Moccasin
                    ]
                },
                {
                    6, [
                        OxyColors.Cornsilk, OxyColors.BlanchedAlmond, OxyColors.Wheat, OxyColors.Goldenrod,
                        OxyColors.DarkKhaki, OxyColors.Khaki, OxyColors.Moccasin
                    ]
                },
                {
                    7, [
                        OxyColors.Cornsilk, OxyColors.BlanchedAlmond, OxyColors.Wheat, OxyColors.Goldenrod,
                        OxyColors.DarkKhaki, OxyColors.Khaki, OxyColors.Moccasin
                    ]
                },
                {
                    8, [
                        OxyColors.Cornsilk, OxyColors.BlanchedAlmond, OxyColors.Wheat, OxyColors.Goldenrod,
                        OxyColors.DarkKhaki, OxyColors.Khaki, OxyColors.Moccasin
                    ]
                },
                {
                    9, [
                        OxyColors.DarkBlue, OxyColors.SkyBlue, OxyColors.CornflowerBlue,
                        OxyColors.DarkBlue, OxyColors.CadetBlue, OxyColors.SteelBlue, OxyColors.DodgerBlue,
                        OxyColors.AliceBlue, OxyColors.DarkSlateBlue
                    ]
                },
                {
                    10, [
                        OxyColors.DarkRed, OxyColors.LightCoral, OxyColors.PaleVioletRed,
                        OxyColors.IndianRed, OxyColors.Firebrick, OxyColors.Maroon, OxyColors.Tomato
                    ]
                },
                {
                    11, [
                        OxyColors.Green, OxyColors.MediumSpringGreen, OxyColors.LightGreen,
                        OxyColors.Linen, OxyColors.SpringGreen, OxyColors.Chartreuse, OxyColors.DarkSeaGreen
                    ]
                }
            };

            IEnumerable<OxyColor> overflow = new List<OxyColor>()
            {
                OxyColors.Cornsilk, OxyColors.BlanchedAlmond, OxyColors.Aqua, OxyColors.Aquamarine, 
                OxyColors.HotPink, OxyColors.PaleGreen, OxyColors.Gray, OxyColors.SeaGreen,
                OxyColors.LemonChiffon, OxyColors.RosyBrown, OxyColors.MediumSpringGreen
            };
            overflowColors = new Queue<OxyColor>(overflow);
        }

    }
}
