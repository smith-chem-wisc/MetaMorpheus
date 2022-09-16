using Chemistry;
using EngineLayer;
using MassSpectrometry;
using OxyPlot;
using Proteomics.Fragmentation;
using System.Collections.Generic;
using System.Linq;
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

        public List<PsmFromTsv> SpectrumMatches { get; private set; }
        public Dictionary<string, List<PsmFromTsv>> PsmsByProteinDictionary { get; private set; }

        public ChimeraSpectrumMatchPlot(PlotView plotView, MsDataScan scan, List<PsmFromTsv> psms) : base(plotView, null, scan)
        {
            SpectrumMatches = psms;
            PsmsByProteinDictionary = SpectrumMatches.GroupBy(p => p.BaseSeq).ToDictionary(p => p.Key, p => p.ToList());
            psms.Select(p => p.MatchedIons).ForEach(p => matchedFragmentIons.AddRange(p));
            
            AnnotateMatchedIons();
            ZoomAxes();
            RefreshChart();
        }

        /// <summary>
        /// Annotates the matched ions based upon the protein of origin, and the unique protoeform ID's
        /// </summary>
        protected new void AnnotateMatchedIons()
        {
            List<MatchedFragmentIon> allMatchedIons = new();
            List<(string, MatchedFragmentIon)> allDrawnIons = new();
            Queue<OxyColor> overflowColors = OverflowColors;

            int proteinIndex = 0;
            foreach (var proteinGroup in PsmsByProteinDictionary.Values)
            {
                List<MatchedFragmentIon> proteinMatchedIons = new();
                List<MatchedFragmentIon> proteinDrawnIons = new();

                for (int j = 0; j < proteinGroup.Count; j++)
                {
                    proteinMatchedIons.AddRange(proteinGroup[j].MatchedIons);
                    allMatchedIons.AddRange(proteinGroup[j].MatchedIons);
                    PeptideWithSetModifications pepWithSetMods = new(proteinGroup[j].FullSequence.Split('|')[0], GlobalVariables.AllModsKnownDictionary);

                    // more proteins than protein programmed colors
                    if (proteinIndex >= ColorByProteinDictionary.Keys.Count)
                    {
                        proteinIndex = 0;
                    }

                    // each matched ion
                    foreach (var matchedIon in proteinGroup[j].MatchedIons)
                    {
                        OxyColor color = MultipleProteinSharedColor;
                        // if drawn by the same protein already
                        if (proteinDrawnIons.Any(p => p.Annotation == matchedIon.Annotation && p.Mz == matchedIon.Mz))
                        {
                            color = ColorByProteinDictionary[proteinIndex][0];
                        }
                        // if unique peak
                        else
                        {
                            // more proteoforms than programmed colors
                            if (j + 1 >= ColorByProteinDictionary[proteinIndex].Count)
                            {
                                color = overflowColors.Dequeue();
                            }
                            else
                            {
                                color = ColorByProteinDictionary[proteinIndex][j + 1];
                            }
                            proteinDrawnIons.Add(matchedIon);
                        }
                        AnnotatePeak(matchedIon, false, false, color);
                        allDrawnIons.Add((proteinGroup[j].BaseSeq, matchedIon));
                    }
                }
                proteinIndex++;
            }
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
                OxyColors.Blue, OxyColors.Navy, OxyColors.SkyBlue, OxyColors.CornflowerBlue,
                OxyColors.DarkBlue, OxyColors.CadetBlue, OxyColors.SteelBlue, OxyColors.DodgerBlue
            });
            ColorByProteinDictionary.Add(1, new List<OxyColor>()
            {
                OxyColors.Red, OxyColors.DarkRed, OxyColors.LightCoral, OxyColors.PaleVioletRed,
                OxyColors.IndianRed, OxyColors.Firebrick, OxyColors.Maroon, OxyColors.Tomato
            });
            ColorByProteinDictionary.Add(2, new List<OxyColor>()
            {
                OxyColors.Green, OxyColors.DarkGreen, OxyColors.MediumSpringGreen, OxyColors.LightGreen,
                OxyColors.Linen, OxyColors.SpringGreen, OxyColors.Chartreuse, OxyColors.DarkSeaGreen
            });
            ColorByProteinDictionary.Add(3, new List<OxyColor>()
            {
                OxyColors.Purple, OxyColors.Indigo, OxyColors.MediumPurple, OxyColors.Violet,
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
