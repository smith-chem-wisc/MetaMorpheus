using Chemistry;
using EngineLayer;
using iText.IO.Image;
using iText.Kernel.Pdf;
using iText.Layout;
using MassSpectrometry;
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
using System.Threading;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using Easy.Common.Extensions;
using iText.Kernel.Geom;
using OxyPlot.Wpf;
using Proteomics.ProteolyticDigestion;
using ThermoFisher.CommonCore.Data.Business;
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
        public static OxyColor MultipleProteinSharedColor;
        public static Dictionary<int, List<OxyColor>> ColorByProteinDictionary;
        public List<PsmFromTsv> SpectrumMatches { get; protected set; }
        public Dictionary<string, List<PsmFromTsv>> PsmsByProteinDictionary { get; protected set; }

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

                    // each matched ion
                    foreach (var matchedIon in proteinGroup[j].MatchedIons)
                    {
                        // if drawn by the same protein already
                        if (proteinDrawnIons.Any(p => p.Annotation == matchedIon.Annotation && p.Mz == matchedIon.Mz))
                        {
                            AnnotatePeak(matchedIon, false, false, ColorByProteinDictionary[proteinIndex][0]);
                        }
                        // if drawn by a different protein
                        else if (allDrawnIons.Any(p => p.Item1 != proteinGroup[j].BaseSeq && p.Item2.Annotation == matchedIon.Annotation && p.Item2.Mz == matchedIon.Mz))
                        {
                            AnnotatePeak(matchedIon, false, false, MultipleProteinSharedColor);
                        }
                        // if unique peak
                        else
                        {
                            AnnotatePeak(matchedIon, false, false, ColorByProteinDictionary[proteinIndex][j + 1]);
                            proteinDrawnIons.Add(matchedIon);
                        }
                        allDrawnIons.Add((proteinGroup[j].BaseSeq, matchedIon));
                    }
                }
                proteinIndex++;
            }

            //// each protein
            //for (int i = 0; i < PsmsByProteinDictionary.Count; i++)
            //{
            //    var matchedIonsLists = PsmsByProteinDictionary[i].Select(p => p.MatchedIons).ToList();
            //    List<MatchedFragmentIon> proteinMatchedIons = new();
            //    List<MatchedFragmentIon> proteinDrawnIons = new();

            //    // each proteoform
            //    int proteoformIndex = 1;
            //    foreach (var proteoform in proteinGroups[i])
            //    {
            //        proteinMatchedIons.AddRange(proteoform.MatchedIons);
            //        allMatchedIons.AddRange(proteoform.MatchedIons);
            //        PeptideWithSetModifications pepWithSetMods = new(proteoform.FullSequence.Split('|')[0], GlobalVariables.AllModsKnownDictionary);

            //        // each matched ion
            //        foreach (var matchedIon in proteoform.MatchedIons)
            //        {
            //            // if drawn by the same protein alread
            //            if (proteinDrawnIons.Any(p => p.Annotation == matchedIon.Annotation && p.Mz == matchedIon.Mz))
            //                AnnotatePeak(matchedIon, false, false, SingleProteinColors.ColorsByProteinDict[i][0].Item1);
            //            // if drawn by a different protein
            //            else if (allDrawnIons.Any(p => p.Item1 != proteoform.BaseSeq && p.Item2.Annotation == matchedIon.Annotation && p.Item2.Mz == matchedIon.Mz))

            //                AnnotatePeak(matchedIon, false, false, OxyColors.Black);
            //            // if unique peak
            //            else
            //                AnnotatePeak(matchedIon, false, false, SingleProteinColors.ColorsByProteinDict[i][proteoformIndex].Item1);
            //            proteinDrawnIons.Add(matchedIon);
            //            allDrawnIons.Add((proteoform.BaseSeq, matchedIon));
            //        }
            //        proteoformIndex++;
            //    }
            //}
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
                OxyColors.Red, OxyColors.DarkRed, OxyColors.LightCoral, OxyColors.PaleVioletRed,
                OxyColors.IndianRed, OxyColors.Firebrick, OxyColors.Maroon, OxyColors.Tomato
            });
            ColorByProteinDictionary.Add(1, new List<OxyColor>()
            {
                OxyColors.Blue, OxyColors.Navy, OxyColors.SkyBlue, OxyColors.CornflowerBlue,
                OxyColors.DarkBlue, OxyColors.CadetBlue, OxyColors.SteelBlue, OxyColors.DodgerBlue
            });
            ColorByProteinDictionary.Add(2, new List<OxyColor>()
            {
                OxyColors.Purple, OxyColors.Indigo, OxyColors.MediumPurple, OxyColors.Violet,
                OxyColors.Plum, OxyColors.Orchid, OxyColors.BlueViolet, OxyColors.Magenta
            });
            ColorByProteinDictionary.Add(3, new List<OxyColor>()
            {
                OxyColors.Green, OxyColors.DarkGreen, OxyColors.MediumSpringGreen, OxyColors.LightGreen,
                OxyColors.Linen, OxyColors.SpringGreen, OxyColors.Chartreuse, OxyColors.DarkSeaGreen
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
        }

    }
}
