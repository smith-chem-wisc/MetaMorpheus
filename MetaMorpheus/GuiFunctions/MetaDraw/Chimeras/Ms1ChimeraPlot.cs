using MzLibUtil;
using OxyPlot.Wpf;

namespace GuiFunctions.MetaDraw;

public class Ms1ChimeraPlot : SpectrumMatchPlot
{
    public ChimeraGroupViewModel ChimeraGroup { get; }
    public DoubleRange Range { get; }

    public Ms1ChimeraPlot(PlotView plotView, ChimeraGroupViewModel chimeraGroupVm) : base(plotView, null,
        chimeraGroupVm.Ms1Scan)
    {
        Range = chimeraGroupVm.Ms2Scan.IsolationRange;
        ChimeraGroup = chimeraGroupVm;

        ZoomAxes();
        AnnotateChimericPeaks(chimeraGroupVm);
        AnnotateIsolationWindow(Range);
        RefreshChart();
    }

    private void AnnotateChimericPeaks(ChimeraGroupViewModel chimeraGroupVm)
    {
        foreach (var ionGroup in chimeraGroupVm.PrecursorIonsByColor)
        {
            var color = ionGroup.Key;
            var ions = ionGroup.Value;
            foreach (var tuple in ions)
            {
                AnnotatePeak(tuple.Item1, false, false, color, tuple.Item2);
            }
        }
    }

    public void ZoomAxes()
    {
        Model.Axes[0].MajorStep = 1;
        Model.Axes[0].MinorStep = 0.2;

        var extracted = ChimeraGroup.Ms1Scan.MassSpectrum.Extract(Range);
        double maxIntensity = 0;
        foreach (var p in extracted)
        {
            if (p.Intensity > maxIntensity)
                maxIntensity = p.Intensity;
        }
        maxIntensity *= 1.4;

        Model.Axes[0].Zoom(Range.Minimum - 1, Range.Maximum + 1);
        Model.Axes[1].Zoom(0, maxIntensity);
    }
}