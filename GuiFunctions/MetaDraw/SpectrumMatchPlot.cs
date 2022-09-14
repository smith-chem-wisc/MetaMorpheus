using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using EngineLayer;
using MassSpectrometry;
using mzPlot;
using OxyPlot;
using OxyPlot.Annotations;
using OxyPlot.Axes;
using OxyPlot.Series;
using Proteomics.Fragmentation;

namespace GuiFunctions
{
    public class SpectrumMatchPlot : Plot
    {
        private List<MatchedFragmentIon> matchedFragmentIons;
        public MsDataScan Scan { get; protected set; }
        public PsmFromTsv SpectrumMatch { get; set; }

        public SpectrumMatchPlot(OxyPlot.Wpf.PlotView plotView, PsmFromTsv psm, MsDataScan scan) : base(plotView)
        {
            Model.Title = string.Empty;
            Model.Subtitle = string.Empty;
            Scan = scan;
            DrawSpectrum();

            if (psm != null)
            {
                SpectrumMatch = psm;
                matchedFragmentIons = SpectrumMatch.MatchedIons;
                AnnotateMatchedIons(isBetaPeptide: false, matchedFragmentIons);
            }
            RefreshChart();
        }

        /// <summary>
        /// 
        /// </summary>
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
                AxisTitleDistance = 10,
                ExtraGridlines = new double[] { 0 },
                ExtraGridlineColor = OxyColors.Black,
                ExtraGridlineThickness = 1
            });

            // draw all peaks in the scan
            for (int i = 0; i < Scan.MassSpectrum.XArray.Length; i++)
            {
                double mz = Scan.MassSpectrum.XArray[i];
                double intensity = Scan.MassSpectrum.YArray[i];

                DrawPeak(mz, intensity, MetaDrawSettings.StrokeThicknessUnannotated, MetaDrawSettings.UnannotatedPeakColor, null);
            }
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

        protected void AnnotateMatchedIons(bool isBetaPeptide, List<MatchedFragmentIon> matchedFragmentIons, bool useLiteralPassedValues = false)
        {
            foreach (MatchedFragmentIon matchedIon in matchedFragmentIons)
            {
                AnnotatePeak(matchedIon, isBetaPeptide, useLiteralPassedValues);
            }
        }

        protected void AnnotatePeak(MatchedFragmentIon matchedIon, bool isBetaPeptide, bool useLiteralPassedValues = false)
        {
            OxyColor ionColor;

            if (SpectrumMatch.VariantCrossingIons.Contains(matchedIon))
            {
                ionColor = MetaDrawSettings.VariantCrossColor;
            }
            else if (matchedIon.NeutralTheoreticalProduct.SecondaryProductType != null) //if internal fragment
            {
                ionColor = MetaDrawSettings.InternalIonColor;
            }
            else if (isBetaPeptide)
            {
                ionColor = MetaDrawSettings.BetaProductTypeToColor[matchedIon.NeutralTheoreticalProduct.ProductType];
            }
            else
            {
                ionColor = MetaDrawSettings.ProductTypeToColor[matchedIon.NeutralTheoreticalProduct.ProductType];
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
}
