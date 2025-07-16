using System.Collections.Generic;
using System.Linq;
using MzLibUtil;
using OxyPlot;
using OxyPlot.Wpf;
using LineSeries = OxyPlot.Series.LineSeries;

namespace GuiFunctions.MetaDraw.Chimeras;

public class Ms1ChimeraPlot : SpectrumMatchPlot
{
    public ChimeraGroupViewModel ChimeraGroup { get; private set; }
    public PlotView PlotView { get; private set; }
    public DoubleRange Range { get; private set; }

    public Ms1ChimeraPlot(PlotView plotView, ChimeraGroupViewModel chimeraGroupVm) : base(plotView, null,
        chimeraGroupVm.Ms1Scan)
    {
        PlotView = plotView;
        Range = chimeraGroupVm.Ms2Scan.IsolationRange;
        ChimeraGroup = chimeraGroupVm;

        ZoomAxes();
        AnnotateChimericPeaks(chimeraGroupVm);
        SetTitle();
        AnnotateIsolationWindow();
        RefreshChart();
    }

    private void AnnotateIsolationWindow()
    {
        var isolationWindow = ChimeraGroup.Ms2Scan.IsolationRange;
        List<DataPoint> points = new List<DataPoint>()
            {
                new(isolationWindow.Minimum, Model.Axes[1].AbsoluteMaximum),
                new(isolationWindow.Minimum, 0),
                new(isolationWindow.Maximum, 0),
                new(isolationWindow.Maximum, Model.Axes[1].AbsoluteMaximum),
            };
        var lineSeries = new LineSeries();
        lineSeries.Points.AddRange(points);
        lineSeries.LineStyle = LineStyle.Dash;
        lineSeries.Color = OxyColors.Red;
        Model.Series.Add(lineSeries);
    }

    private void AnnotateChimericPeaks(ChimeraGroupViewModel chimeraGroupVm)
    {
        foreach (var ionGroup in chimeraGroupVm.PrecursorIonsByColor)
        {
            var color = ionGroup.Key;
            ionGroup.Value.ForEach(p => AnnotatePeak(p.Item1, false, false, color, p.Item2));
        }
    }

    private void SetTitle()
    {
        string title = "";


        Model.Title = title;
    }

    public void ZoomAxes()
    {
        Model.Axes[0].MajorStep = 1;
        Model.Axes[0].MinorStep = 0.2;
        var maxIntensity = ChimeraGroup.Ms1Scan.MassSpectrum.Extract(Range).Max(p => p.Intensity) * 1.4;
        Model.Axes[0].Zoom(Range.Minimum - 1, Range.Maximum + 1);
        Model.Axes[1].Zoom(0, maxIntensity);
    }
}