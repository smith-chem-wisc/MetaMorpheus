using Chemistry;
using EngineLayer;
using MassSpectrometry;
using OxyPlot;
using Omics.Fragmentation;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using Easy.Common.Extensions;
using OxyPlot.Wpf;
using Proteomics.ProteolyticDigestion;
using Point = System.Windows.Point;
using Vector = System.Windows.Vector;
using Canvas = System.Windows.Controls.Canvas;
using LinearAxis = OxyPlot.Axes.LinearAxis;
using LineSeries = OxyPlot.Series.LineSeries;
using Plot = mzPlot.Plot;
using TextAnnotation = OxyPlot.Annotations.TextAnnotation;
using Readers;

namespace GuiFunctions
{
    public class ChimeraSpectrumMatchPlot : SpectrumMatchPlot
    {
        private static Queue<OxyColor> overflowColors;
        public static OxyColor MultipleProteinSharedColor;
        public static Dictionary<int, List<OxyColor>> ColorByProteinDictionary;
        public static Queue<OxyColor> OverflowColors
        {
            get => new Queue<OxyColor>(overflowColors.ToList());
        }

        public List<SpectrumMatchFromTsv> SpectrumMatches { get; private set; }
        public Dictionary<string, List<SpectrumMatchFromTsv>> PsmsByProteinDictionary { get; private set; }

        public ChimeraSpectrumMatchPlot(PlotView plotView, MsDataScan scan, List<SpectrumMatchFromTsv> sms) : base(plotView, null, scan)
        {
            SpectrumMatches = sms;
            PsmsByProteinDictionary = SpectrumMatches.GroupBy(p => p.BaseSeq).ToDictionary(p => p.Key, p => p.ToList());
            sms.Select(p => p.MatchedIons).ForEach(p => MatchedFragmentIons.AddRange(p));
            
            AnnotateMatchedIons();
            ZoomAxes();
            RefreshChart();
        }

        /// <summary>
        /// Annotates the matched ions based upon the protein of origin, and the unique proteoform ID's
        /// </summary>
        protected void AnnotateMatchedIons()
        {
            var allMatchedIons = new List<MatchedFragmentIon>();
            List<(string, MatchedFragmentIon)> allDrawnIons = new();
            var overflowColors = OverflowColors;

            var proteinIndex = 0;
            foreach (var proteinGroup in PsmsByProteinDictionary.Values)
            {
                var proteinMatchedIons = new List<MatchedFragmentIon>();
                var proteinDrawnIons = new List<MatchedFragmentIon>();

                for (var j = 0; j < proteinGroup.Count; j++)
                {
                    proteinMatchedIons.AddRange(proteinGroup[j].MatchedIons);
                    allMatchedIons.AddRange(proteinGroup[j].MatchedIons);
                    var bioPolymerWithSetMods = proteinGroup.First()
                        .ToBioPolymerWithSetMods(proteinGroup[j].FullSequence.Split('|')[0]);

                    // more proteins than protein programmed colors
                    if (proteinIndex >= ColorByProteinDictionary.Keys.Count)
                    {
                        proteinIndex = 0;
                    }

                    // each matched ion
                    foreach (var matchedIon in proteinGroup[j].MatchedIons)
                    {
                        if (!MetaDrawSettings.DisplayInternalIons &&
                            matchedIon.NeutralTheoreticalProduct.SecondaryProductType != null)
                            continue;

                        OxyColor color;

                        // if drawn by the same protein already
                        if (proteinDrawnIons.Any(p => p.Equals(matchedIon)))
                        {
                            color = ColorByProteinDictionary[proteinIndex][0];
                        }
                        // if drawn already by different protein
                        else if (allDrawnIons.Any(p => p.Item2.Equals(matchedIon)))
                        {
                            color = MultipleProteinSharedColor;
                        }
                        // if unique peak
                        else
                        {
                            // more proteoforms than programmed colors
                            if (j + 1 >= ColorByProteinDictionary[proteinIndex].Count)
                            {
                                ColorByProteinDictionary[proteinIndex].Add(overflowColors.Dequeue());
                            }
                            color = ColorByProteinDictionary[proteinIndex][j + 1];
                            proteinDrawnIons.Add(matchedIon);
                        }
                        AnnotatePeak(matchedIon, false, false, color);
                        allDrawnIons.Add((proteinGroup[j].BaseSeq, matchedIon));
                    }
                }
                proteinIndex++;
            }
        }

        public void ExportPlot(string path, Canvas legend = null,  double width = 700, double height = 370)
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

            // render legend as bitmap and export as png if used
            System.Drawing.Bitmap ptmLegendBitmap = null;
            Point legendPoint;
            if (legend != null && MetaDrawSettings.ShowLegend)
            {
                // Saving Canvas as a usable Png
                RenderTargetBitmap legendRenderBitmap = new((int)(dpiScale * legend.ActualWidth), (int)(dpiScale * legend.ActualHeight),
                         MetaDrawSettings.CanvasPdfExportDpi, MetaDrawSettings.CanvasPdfExportDpi, PixelFormats.Pbgra32);
                legendRenderBitmap.Render(legend);
                var legendEncoder = new PngBitmapEncoder();
                legendEncoder.Frames.Add(BitmapFrame.Create(legendRenderBitmap));
                using (var file = File.Create(tempLegendPngPath))
                {
                    legendEncoder.Save(file);
                }

                // converting png to the final bitmap format
                System.Drawing.Bitmap tempLegendBitmap = new(tempLegendPngPath);
                ptmLegendBitmap = new System.Drawing.Bitmap(tempLegendBitmap, new System.Drawing.Size((int)legend.ActualWidth, (int)legend.ActualHeight));
                bitmaps.Add(ptmLegendBitmap);
                legendPoint = new Point(0, height);
                points.Add(legendPoint);
                base.ExportPlot(path, ptmLegendBitmap, width, height);
                tempLegendBitmap.Dispose();
            }

            // combine the bitmaps
            var combinedBitmaps = MetaDrawLogic.CombineBitmap(bitmaps, points, false);
            File.Delete(tempModelPath);
            File.Delete(tempLegendPngPath);
            base.ExportPlot(path, combinedBitmaps, width, height);
        }

        /// <summary>
        /// Initializes the colors to be used by the Chimera Plotter
        /// </summary>
        static ChimeraSpectrumMatchPlot()
        {
            MultipleProteinSharedColor = OxyColors.Black;
            ColorByProteinDictionary = new();
            ColorByProteinDictionary.Add(0, new List<OxyColor>()
            {
                OxyColors.Blue, OxyColors.SkyBlue, OxyColors.CornflowerBlue,
                OxyColors.DarkBlue, OxyColors.CadetBlue, OxyColors.SteelBlue, OxyColors.DodgerBlue
            });
            ColorByProteinDictionary.Add(1, new List<OxyColor>()
            {
                OxyColors.Red, OxyColors.LightCoral, OxyColors.PaleVioletRed,
                OxyColors.IndianRed, OxyColors.Firebrick, OxyColors.Maroon, OxyColors.Tomato
            });
            ColorByProteinDictionary.Add(2, new List<OxyColor>()
            {
                OxyColors.Green, OxyColors.MediumSpringGreen, OxyColors.LightGreen,
                OxyColors.Linen, OxyColors.SpringGreen, OxyColors.Chartreuse, OxyColors.DarkSeaGreen
            });
            ColorByProteinDictionary.Add(3, new List<OxyColor>()
            {
                OxyColors.Purple, OxyColors.MediumPurple, OxyColors.Violet,
                OxyColors.Plum, OxyColors.Orchid, OxyColors.BlueViolet, OxyColors.Magenta
            });
            ColorByProteinDictionary.Add(4, new List<OxyColor>()
            {
                OxyColors.Brown, OxyColors.SaddleBrown, OxyColors.Sienna, OxyColors.Chocolate,
                OxyColors.SandyBrown, OxyColors.Chocolate, OxyColors.Peru, OxyColors.Tan
            });
            ColorByProteinDictionary.Add(5, new List<OxyColor>()
            {
                OxyColors.Gold, OxyColors.DarkGoldenrod, OxyColors.Wheat, OxyColors.Goldenrod,
                OxyColors.DarkKhaki, OxyColors.Khaki, OxyColors.Moccasin
            });

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
