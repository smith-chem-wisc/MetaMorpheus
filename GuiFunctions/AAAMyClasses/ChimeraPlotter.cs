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
    public class ChimeraPlotter : PeptideSpectrumMatchPlot
    {
        public List<PsmFromTsv> SpectrumMatches { get; protected set; }
        public List<PtmLegendViewModel> Legends { get; protected set; }

        public ChimeraPlotter(PlotView plotView, MsDataScan scan, List<PsmFromTsv> psms) : base(plotView, null, scan, null, false)
        {
            Legends = new();
            SpectrumMatches = psms;
            Scan = scan;

            List<MatchedFragmentIon> allMatchedIons = new();
            List<(string, MatchedFragmentIon)> allDrawnIons = new();
            var proteinGroups = SpectrumMatches.GroupBy(p => p.BaseSeq).ToList();

            // each protein
            for (int i = 0; i < proteinGroups.Count; i++)
            {
                var matchedIonsLists = proteinGroups[i].Select(p => p.MatchedIons).ToList();
                List<MatchedFragmentIon> proteinMatchedIons = new();
                List<MatchedFragmentIon> proteinDrawnIons = new();
                List<ChimeraLegendItemViewModel> legendItems = new();
                legendItems.Add(new("Shared Ions", ChimeraPlotterColors.ColorsByProteinDict[i][0]));

                // each proteoform
                int proteoformIndex = 1;
                foreach (var proteoform in proteinGroups[i])
                {
                    proteinMatchedIons.AddRange(proteoform.MatchedIons);
                    allMatchedIons.AddRange(proteoform.MatchedIons);
                    PeptideWithSetModifications pepWithSetMods = new(proteoform.FullSequence.Split('|')[0], GlobalVariables.AllModsKnownDictionary);
                    var legendItem = new ChimeraLegendItemViewModel(String.Join(", ",
                            pepWithSetMods.AllModsOneIsNterminus.Select(p => p.Key + " - " + p.Value.IdWithMotif)
                                .ToArray()), ChimeraPlotterColors.ColorsByProteinDict[i][proteoformIndex]);
                    legendItems.Add(legendItem);

                    // each matched ion
                    foreach (var matchedIon in proteoform.MatchedIons)
                    {
                        // if drawn by the same protein alread
                        if (proteinDrawnIons.Any(p => p.Annotation == matchedIon.Annotation && p.Mz == matchedIon.Mz))
                            AnnotatePeak(matchedIon, false, ChimeraPlotterColors.ColorsByProteinDict[i][0].Item1,
                                ChimeraPlotterColors.ColorsByProteinDict[i][0].Item2, ChimeraPlotterColors.ColorsByProteinDict[i][0].Item3);
                        // if drawn by a different protein
                        else if (allDrawnIons.Any(p => p.Item1 != proteoform.BaseSeq && p.Item2.Annotation == matchedIon.Annotation && p.Item2.Mz == matchedIon.Mz))

                            AnnotatePeak(matchedIon, false, OxyColors.Black, OxyColors.Black, OxyColors.Black);
                        // if unique peak
                        else
                            AnnotatePeak(matchedIon, false, ChimeraPlotterColors.ColorsByProteinDict[i][proteoformIndex].Item1,
                                ChimeraPlotterColors.ColorsByProteinDict[i][proteoformIndex].Item2, ChimeraPlotterColors.ColorsByProteinDict[i][proteoformIndex].Item3);
                        proteinDrawnIons.Add(matchedIon);
                        allDrawnIons.Add((proteoform.BaseSeq, matchedIon));
                    }
                    proteoformIndex++;
                }
                Legends.Add(new(legendItems, proteinGroups[i].First().BaseSeq));
            }

           

            ZoomAxes(allMatchedIons);
            RefreshChart();
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


        protected void AnnotatePeak(MatchedFragmentIon matchedIon, bool isBetaPeptide, OxyColor bColor, OxyColor yColor, OxyColor internalColor, bool useLiteralPassedValues = false)
        {
            OxyColor ionColor = OxyColors.Black;

            if (matchedIon.NeutralTheoreticalProduct.SecondaryProductType != null) //if internal fragment
            {
                if (MetaDrawSettings.InternalIonColor == OxyColors.Transparent)
                    ionColor = OxyColors.Transparent;
                else
                    ionColor = bColor;
            }
            else if (matchedIon.NeutralTheoreticalProduct.ProductType == ProductType.b)
            {
                ionColor = bColor;
            }
            else if (matchedIon.NeutralTheoreticalProduct.ProductType == ProductType.y)
            {
                ionColor = bColor;
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
            var peakAnnotation = new TextAnnotation();
            if (MetaDrawSettings.DisplayIonAnnotations)
            {
                string peakAnnotationText = prefix + matchedIon.NeutralTheoreticalProduct.Annotation;

                if (matchedIon.NeutralTheoreticalProduct.NeutralLoss != 0 && !peakAnnotationText.Contains("-" + matchedIon.NeutralTheoreticalProduct.NeutralLoss.ToString("F2")))
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


                peakAnnotation.Font = "Arial";
                peakAnnotation.FontSize = MetaDrawSettings.AnnotatedFontSize;
                peakAnnotation.FontWeight = MetaDrawSettings.AnnotationBold ? OxyPlot.FontWeights.Bold : 2.0;
                peakAnnotation.TextColor = ionColor;
                peakAnnotation.StrokeThickness = 0;
                peakAnnotation.Text = peakAnnotationText;
                peakAnnotation.TextPosition = new DataPoint(mz, intensity);
                peakAnnotation.TextVerticalAlignment = intensity < 0 ? OxyPlot.VerticalAlignment.Top : OxyPlot.VerticalAlignment.Bottom;
                peakAnnotation.TextHorizontalAlignment = OxyPlot.HorizontalAlignment.Center;
            }
            else
            {
                peakAnnotation.Text = string.Empty;
            }
            if (matchedIon.NeutralTheoreticalProduct.SecondaryProductType != null && !MetaDrawSettings.DisplayInternalIonAnnotations) //if internal fragment
            {
                peakAnnotation.Text = string.Empty;
            }

            DrawPeak(mz, intensity, MetaDrawSettings.StrokeThicknessAnnotated, ionColor, peakAnnotation);
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

    public static class ChimeraPlotterColors
    {
        public static Dictionary<int, List<(OxyColor, OxyColor, OxyColor)>> ColorsByProteinDict { get; set; }
        public static OxyColor SharedBetweenProteinColor { get; set; }

        static ChimeraPlotterColors()
        {
            InitializeColorValues();
        }

        private static void InitializeColorValues()
        {
            SharedBetweenProteinColor = OxyColors.Black;
            List<(OxyColor, OxyColor, OxyColor)> proteinAColors = new List<(OxyColor, OxyColor, OxyColor)>();
            List<(OxyColor, OxyColor, OxyColor)> proteinBColors = new List<(OxyColor, OxyColor, OxyColor)>();
            List<(OxyColor, OxyColor, OxyColor)> proteinCColors = new List<(OxyColor, OxyColor, OxyColor)>();
            List<(OxyColor, OxyColor, OxyColor)> proteinDColors = new List<(OxyColor, OxyColor, OxyColor)>();


            // first color is that of shared peak for the protein, rest are for unique peaks to the proteoforms
            proteinBColors.Add((OxyColors.DarkBlue, OxyColors.DarkBlue, OxyColors.LightBlue));
            proteinBColors.Add((OxyColors.LightCyan, OxyColors.DarkCyan, OxyColors.LightCyan));
            proteinBColors.Add((OxyColors.LightSkyBlue, OxyColors.DeepSkyBlue, OxyColors.LightSkyBlue));
            proteinBColors.Add((OxyColors.LightSteelBlue, OxyColors.RoyalBlue, OxyColors.LightSteelBlue));
            proteinBColors.Add((OxyColors.PaleTurquoise, OxyColors.DarkTurquoise, OxyColors.PaleTurquoise));

            proteinCColors.Add((OxyColors.DarkGreen, OxyColors.DarkGreen, OxyColors.LightGreen));
            proteinCColors.Add((OxyColors.LightSeaGreen, OxyColors.DarkSeaGreen, OxyColors.LightSeaGreen));
            proteinCColors.Add((OxyColors.LimeGreen, OxyColors.Chartreuse, OxyColors.LimeGreen));
            proteinCColors.Add((OxyColors.Aquamarine, OxyColors.MediumAquamarine, OxyColors.Aquamarine));
            proteinCColors.Add((OxyColors.YellowGreen, OxyColors.LawnGreen, OxyColors.YellowGreen));

            proteinDColors.Add((OxyColors.Indigo, OxyColors.Indigo, OxyColors.Violet));
            proteinDColors.Add((OxyColors.Red, OxyColors.DeepPink, OxyColors.LightPink));
            proteinDColors.Add((OxyColors.Blue, OxyColors.Fuchsia, OxyColors.PaleVioletRed));
            proteinDColors.Add((OxyColors.MediumOrchid, OxyColors.DarkOrchid, OxyColors.MediumOrchid));
            proteinDColors.Add((OxyColors.MediumPurple, OxyColors.DarkMagenta, OxyColors.MediumPurple));

            proteinAColors.Add((OxyColors.SaddleBrown, OxyColors.SaddleBrown, OxyColors.BurlyWood));
            proteinAColors.Add((OxyColors.LightCoral, OxyColors.DarkRed, OxyColors.LightCoral));
            proteinAColors.Add((OxyColors.Peru, OxyColors.DarkOrange, OxyColors.Peru));
            proteinAColors.Add((OxyColors.LightSalmon, OxyColors.DarkSalmon, OxyColors.LightSalmon));
            proteinAColors.Add((OxyColors.PaleVioletRed, OxyColors.Fuchsia, OxyColors.PaleVioletRed));

            ColorsByProteinDict = new Dictionary<int, List<(OxyColor, OxyColor, OxyColor)>>();
            ColorsByProteinDict.Add(0, proteinDColors);
            ColorsByProteinDict.Add(1, proteinCColors);
            ColorsByProteinDict.Add(2, proteinBColors);
            ColorsByProteinDict.Add(3, proteinAColors);
        }


    }
}
