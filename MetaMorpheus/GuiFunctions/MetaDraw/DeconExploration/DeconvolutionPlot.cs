#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using GuiFunctions.Util;
using MassSpectrometry;
using MzLibUtil;
using OxyPlot;
using OxyPlot.Wpf;
using FontWeights = OxyPlot.FontWeights;
using HorizontalAlignment = OxyPlot.HorizontalAlignment;
using TextAnnotation = OxyPlot.Annotations.TextAnnotation;
using VerticalAlignment = OxyPlot.VerticalAlignment;

namespace GuiFunctions;

public class DeconvolutionPlot : MassSpectrumPlot
{
    public static readonly CyclicalQueue<OxyColor> ColorQueue = new CyclicalQueue<OxyColor>(new[]
    {
        OxyColors.Red, OxyColors.Blue, OxyColors.Green, OxyColors.Orange, OxyColors.Purple,
        OxyColors.Teal, OxyColors.Brown, OxyColors.Pink, OxyColors.Yellow, OxyColors.Gray,
        OxyColors.Cyan, OxyColors.Magenta, OxyColors.LimeGreen, OxyColors.DarkBlue, OxyColors.DarkRed,
        OxyColors.DarkGreen, OxyColors.Gold, OxyColors.Indigo, OxyColors.Olive, OxyColors.Maroon,
        OxyColors.Navy, OxyColors.Turquoise, OxyColors.Violet, OxyColors.Sienna, OxyColors.Salmon,
        OxyColors.Coral, OxyColors.Khaki, OxyColors.Plum, OxyColors.Peru, OxyColors.SteelBlue,
        OxyColors.MediumPurple, OxyColors.MediumSeaGreen, OxyColors.MediumSlateBlue, OxyColors.MediumVioletRed,
        OxyColors.MediumOrchid, OxyColors.MediumTurquoise, OxyColors.MediumSpringGreen, OxyColors.MediumAquamarine
    });

    public DeconvolutionPlot(PlotView plotView, MsDataScan scan, List<DeconvolutedSpeciesViewModel> deconResults, DeconvolutionMode mode, MzRange? isolationRange = null) : base(plotView, scan)
    {
        double maxAnnotated = deconResults.Any()
            ? deconResults.Max(p => p.Envelope.Peaks.Max(m => m.mz))
            : isolationRange?.Maximum + 1 ?? scan.MassSpectrum.Range.Maximum;

        double maxIntensity;
        if (isolationRange != null)
            maxIntensity = scan.MassSpectrum.Extract(isolationRange).Max(p => p.Intensity);
        else
            maxIntensity = scan.MassSpectrum.YofPeakWithHighestY ?? 0;

        ZoomAxes(scan, maxAnnotated, maxIntensity, mode, isolationRange);

        AnnotatePlot(deconResults);

        if (isolationRange is not null && mode == DeconvolutionMode.IsolationRegion)
            AnnotateIsolationWindow(isolationRange);
        RefreshChart();
    }

    public void AnnotatePlot(List<DeconvolutedSpeciesViewModel> deconResults)
    {
        if (deconResults is null || deconResults.Count == 0)
            return;

        for (var index = 0; index < deconResults.Count; index++)
        {
            var species = deconResults[index];
            species.Color = ColorQueue.Dequeue();
            TextAnnotation annotation = null;
            foreach (var peak in species.Envelope.Peaks)
            {
                int i = Scan.MassSpectrum.GetClosestPeakIndex(peak.mz);
                double mz = Scan.MassSpectrum.XArray[i];
                double intensity = Scan.MassSpectrum.YArray[i];

                if (Math.Abs(peak.mz - species.MostAbundantMz) < 1e-6 && MetaDrawSettings.DisplayIonAnnotations)
                {
                    annotation = new TextAnnotation
                    {
                        Text = species.Annotation,
                        TextColor = species.Color,
                        FontSize = MetaDrawSettings.AnnotatedFontSize,
                        FontWeight = MetaDrawSettings.AnnotationBold ? FontWeights.Bold : 2.0,
                        TextPosition = new DataPoint(mz, intensity),
                        TextVerticalAlignment = intensity < 0 ? VerticalAlignment.Top : VerticalAlignment.Bottom,
                        TextHorizontalAlignment = HorizontalAlignment.Center,
                        StrokeThickness = 0,
                    };
                }
                DrawPeak(mz, intensity, MetaDrawSettings.StrokeThicknessAnnotated, species.Color, annotation);
            }
        }
    }

    private const double YRangeMultiplier = 1.2;
    public void ZoomAxes(MsDataScan scan, double maxAnnotatedMz, double maxAnnotatedIntensity, DeconvolutionMode mode, MzRange? isolationRange = null)
    {
        Model.Title = $"Scan {scan.OneBasedScanNumber} (MS{scan.MsnOrder})";

        // Full Scan
        if (mode is DeconvolutionMode.FullSpectrum)
        {
            var max = isolationRange?.Maximum ?? maxAnnotatedMz + 50;
            var min = isolationRange?.Minimum ?? scan.MassSpectrum.Range.Minimum;

            Model.Axes[0].Zoom(min, max);
            Model.Axes[1].Zoom(0, maxAnnotatedIntensity * YRangeMultiplier);
        }
        // Isolation Region or trimmed display region
        else
        {
            Model.Title = $"Scan {scan.OneBasedScanNumber} (MS{scan.MsnOrder}) - Isolation Window {isolationRange!.Minimum:F2}-{isolationRange.Maximum:F2} m/z";

            Model.Axes[0].Zoom(isolationRange.Minimum - 1, isolationRange.Maximum + 1);
            Model.Axes[1].Zoom(0, maxAnnotatedIntensity * YRangeMultiplier);
        }
    }
}
