using System.Linq;
using System.Threading;
using GuiFunctions;
using GuiFunctions.MetaDraw.Chimeras;
using NUnit.Framework;
using OxyPlot;
using OxyPlot.Wpf;

namespace Test.MetaDraw;

[TestFixture, Apartment(ApartmentState.STA)]
public class Ms1ChimeraPlotTests
{
    [Test]
    public void Ms1ChimeraPlot_Constructor_InitializesProperties()
    {
        // Arrange
        var chimeraGroup = ChimeraGroupViewModelTests.OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        var plotView = new PlotView();

        // Act
        var ms1Plot = new Ms1ChimeraPlot(plotView, chimeraGroup);

        // Assert
        Assert.That(ms1Plot.PlotView, Is.EqualTo(plotView));
        Assert.That(ms1Plot.ChimeraGroup, Is.EqualTo(chimeraGroup));
        Assert.That(ms1Plot.Range, Is.EqualTo(chimeraGroup.Ms2Scan.IsolationRange));
        Assert.That(ms1Plot.Model, Is.Not.Null);
        Assert.That(ms1Plot.Model.Series, Is.Not.Empty, "Model should have at least one series after construction.");
    }

    [Test]
    public void Ms1ChimeraPlot_AnnotateIsolationWindow_AddsDashedRedBox()
    {
        // Arrange
        var chimeraGroup = ChimeraGroupViewModelTests.OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        var plotView = new PlotView();
        var ms1Plot = new Ms1ChimeraPlot(plotView, chimeraGroup);

        // Act
        // Isolation window is annotated in constructor

        // Assert
        var lineSeries = ms1Plot.Model.Series.OfType<OxyPlot.Series.LineSeries>()
            .FirstOrDefault(ls => ls.LineStyle == LineStyle.Dash && ls.Color == OxyColors.Red);
        Assert.That(lineSeries, Is.Not.Null, "A dashed red line series should be present for the isolation window.");
        Assert.That(lineSeries.Points.Count, Is.EqualTo(4), "Isolation window should be a box (4 points).");
        var isolationRange = chimeraGroup.Ms2Scan.IsolationRange;
        Assert.That(lineSeries.Points.Any(p => p.X == isolationRange.Minimum), Is.True);
        Assert.That(lineSeries.Points.Any(p => p.X == isolationRange.Maximum), Is.True);
    }

    [Test]
    public void Ms1ChimeraPlot_AnnotateChimericPeaks_AddsColoredPeaks()
    {
        // Arrange
        var chimeraGroup = ChimeraGroupViewModelTests.TwoProteinsTwoProteoformChimeraGroup.ChimeraGroup;
        var plotView = new PlotView();

        // Act
        var ms1Plot = new Ms1ChimeraPlot(plotView, chimeraGroup);

        // Assert
        // For each precursor color, there should be at least one annotation/peak with that color
        foreach (var color in chimeraGroup.PrecursorIonsByColor.Keys)
        {
            bool found = ms1Plot.Model.Series
                .OfType<OxyPlot.Series.LineSeries>()
                .Any(ls => ls.Color == color);
            Assert.That(found, Is.True, $"Expected to find a peak with color {color} for precursor ion.");
        }
    }

    [Test]
    public void Ms1ChimeraPlot_ZoomAxes_SetsCorrectAxisRanges()
    {
        // Arrange
        var chimeraGroup = ChimeraGroupViewModelTests.OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        var plotView = new PlotView();
        var ms1Plot = new Ms1ChimeraPlot(plotView, chimeraGroup);

        // Act
        ms1Plot.ZoomAxes();

        // Assert
        var xAxis = ms1Plot.Model.Axes[0];
        var yAxis = ms1Plot.Model.Axes[1];
        var range = chimeraGroup.Ms2Scan.IsolationRange;
        Assert.That(xAxis.ActualMinimum, Is.LessThanOrEqualTo(range.Minimum - 1));
        Assert.That(xAxis.ActualMaximum, Is.GreaterThanOrEqualTo(range.Maximum + 1));
        Assert.That(yAxis.ActualMinimum, Is.EqualTo(0).Within(1e-6));
        // yAxis max is set to 1.4 * max intensity in range, so just check it's positive
        Assert.That(yAxis.ActualMaximum, Is.GreaterThan(0));
    }
}