using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using GuiFunctions;
using GuiFunctions.MetaDraw;
using NUnit.Framework;
using OxyPlot;
using OxyPlot.Axes;
using Readers;

namespace Test.MetaDraw;

/// <summary>
/// Phase 3 guard for the histogram renderer swap. OxyPlot 2.1 removed ColumnSeries, and nothing in
/// 2.2 both stacks and draws vertically, so the bars are now StackedColumnSeries rectangles. The
/// failure mode that matters is silent: BarSeries would draw the same numbers rotated 90 degrees,
/// and no assertion on bin values would notice. So these render the model and assert on the
/// screen-space rectangles OxyPlot actually produced, then tie the bar heights back to the bins
/// SpectrumMatchHistograms computed.
/// </summary>
public class StackedColumnSeriesTests
{
    /// <summary>The 14 entries of PlotNames that are histograms; the other two are line plots.</summary>
    private static IEnumerable<string> HistogramPlotNames => PlotModelStat.PlotNames
        .Where(n => n.StartsWith("Histogram of "));

    private static Dictionary<string, ObservableCollection<SpectrumMatchFromTsv>> LoadFixture()
    {
        string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "VariantCrossTest.psmtsv");
        var matches = SpectrumMatchTsvReader.ReadTsv(path, out _);
        return matches.GroupBy(p => p.FileNameWithoutExtension)
            .ToDictionary(g => g.Key, g => new ObservableCollection<SpectrumMatchFromTsv>(g));
    }

    private static PlotModelStat BuildPlot(string plotName, PlotModelStatParameters parameters = null)
    {
        MetaDrawSettings.ResetSettings();
        var bySourceFile = LoadFixture();
        var all = new ObservableCollection<SpectrumMatchFromTsv>(bySourceFile.Values.SelectMany(v => v));
        return new PlotModelStat(plotName, all, bySourceFile, parameters);
    }

    /// <summary>
    /// Forces a full layout and render pass, which is what populates the axis transforms the
    /// geometry assertions read. SvgExporter lives in OxyPlot.Core, so this runs on every platform -
    /// which is the point, since MetaDrawTest cannot.
    /// </summary>
    private static void Render(PlotModel model)
    {
        var exporter = new SvgExporter { Width = 800, Height = 600, IsDocument = true };
        exporter.ExportToString(model);
    }

    private static List<StackedColumnSeries> SeriesOf(PlotModelStat plot) =>
        plot.Model.Series.OfType<StackedColumnSeries>().ToList();

    [Test]
    public void EveryHistogramUsesStackedColumnSeries()
    {
        foreach (var name in HistogramPlotNames)
        {
            var plot = BuildPlot(name);
            Assert.That(plot.Model.Series, Is.Not.Empty, $"{name}: no series");
            Assert.That(plot.Model.Series.All(s => s is StackedColumnSeries), Is.True,
                $"{name}: expected only StackedColumnSeries, got " +
                string.Join(", ", plot.Model.Series.Select(s => s.GetType().Name).Distinct()));
        }
    }

    /// <summary>
    /// The category axis carries the bins and sits on the bottom; the value axis is on the left.
    /// A horizontal bar chart would have these the other way round.
    /// </summary>
    [Test]
    public void CategoryAxisIsOnTheBottomAndValueAxisOnTheLeft()
    {
        foreach (var name in HistogramPlotNames)
        {
            var plot = BuildPlot(name);
            var categoryAxes = plot.Model.Axes.OfType<CategoryAxis>().ToList();
            Assert.That(categoryAxes, Is.Not.Empty, $"{name}: no category axis");
            Assert.That(categoryAxes.All(a => a.Position == AxisPosition.Bottom), Is.True,
                $"{name}: a category axis is not on the bottom");

            var valueAxes = plot.Model.Axes.Where(a => a is not CategoryAxis).ToList();
            Assert.That(valueAxes, Is.Not.Empty, $"{name}: no value axis");
            Assert.That(valueAxes.All(a => a.Position == AxisPosition.Left), Is.True,
                $"{name}: a value axis is not on the left");
        }
    }

    /// <summary>
    /// The rotation check, on real rendered geometry. Bars in different categories must sit at
    /// different horizontal positions and share a common baseline; if the chart were rotated they
    /// would instead stack down the y-axis and share a common left edge.
    /// </summary>
    [Test]
    public void RenderedBarsAreVerticalAndShareABaseline()
    {
        foreach (var name in HistogramPlotNames)
        {
            var plot = BuildPlot(name);
            Render(plot.Model);

            // The axes hold the transforms the render pass just used, so this is the geometry
            // OxyPlot actually drew, not a restatement of the item values.
            var xAxis = plot.Model.Axes.First(a => a.Position == AxisPosition.Bottom);
            var yAxis = plot.Model.Axes.First(a => a.Position == AxisPosition.Left);

            // Bottom-most series only: with stacking, upper segments legitimately float.
            var items = SeriesOf(plot).First().Items.OfType<StackedColumnItem>()
                .Where(i => i.Y1 > i.Y0).ToList();
            if (items.Count < 2)
                continue; // a single-bar histogram cannot show an orientation difference

            var bars = items.Select(i => new
            {
                Left = xAxis.Transform(i.X0),
                Right = xAxis.Transform(i.X1),
                Top = yAxis.Transform(i.Y1),
                Bottom = yAxis.Transform(i.Y0),
            }).ToList();

            Assert.That(bars.Select(b => Round(b.Left)).Distinct().Count(), Is.GreaterThan(1),
                $"{name}: every bar has the same left edge - the chart is rotated");
            Assert.That(bars.Select(b => Round(b.Bottom)).Distinct().Count(), Is.EqualTo(1),
                $"{name}: bars do not share a baseline - the chart is rotated");
            // Screen y grows downward, so a bar's top is above its baseline.
            Assert.That(bars.All(b => b.Top < b.Bottom), Is.True,
                $"{name}: a bar grows downward from its baseline");
            Assert.That(bars.All(b => b.Right > b.Left), Is.True,
                $"{name}: a bar has no horizontal extent");
        }
    }

    /// <summary>
    /// Bars occupy the same fraction of their slot ColumnSeries gave them, i.e. 1 / (1 + GapWidth).
    /// Getting this wrong does not move a number, it just makes every bar the wrong width.
    /// </summary>
    [Test]
    public void BarWidthMatchesTheCategoryAxisGap()
    {
        var plot = BuildPlot("Histogram of Precursor Charges");
        var axis = plot.Model.Axes.OfType<CategoryAxis>().First();
        var series = SeriesOf(plot).First();

        foreach (var item in series.Items.OfType<StackedColumnItem>())
        {
            Assert.That(item.X1 - item.X0, Is.EqualTo(1.0 / (1 + axis.GapWidth)).Within(1e-9));
        }
    }

    /// <summary>
    /// The category axis must span whole slots, -0.5 to N-0.5. A CategoryAxis used to get that from
    /// the categorized series drawn against it; RectangleBarSeries is a plain XYAxisSeries, so left
    /// alone the axis stops at the last bar's right edge and every bar and gridline shifts. Measured
    /// against OxyPlot 2.0: without this the bars come out 4% wide of where ColumnSeries put them.
    /// </summary>
    [Test]
    public void CategoryAxisSpansWholeSlots()
    {
        foreach (var name in HistogramPlotNames)
        {
            var plot = BuildPlot(name);
            Render(plot.Model);

            var axis = plot.Model.Axes.OfType<CategoryAxis>().First(a => a.Key != "GroupAxis");
            int slots = axis.ItemsSource.Cast<object>().Count();

            Assert.That(axis.ActualMinimum, Is.EqualTo(-0.5).Within(1e-9),
                $"{name}: category axis does not start at the left edge of the first slot");
            Assert.That(axis.ActualMaximum, Is.EqualTo(slots - 0.5).Within(1e-9),
                $"{name}: category axis stops at the last bar rather than the end of its slot");
        }
    }

    /// <summary>
    /// Segments of one category stack: each starts where the one below it ended, and the first
    /// starts at the series base value. This is what IsStacked did before it had to be explicit.
    /// </summary>
    [Test]
    public void SegmentsInACategoryAreContiguousFromTheBaseValue()
    {
        var plot = BuildPlot("Histogram of Precursor Charges");
        var series = SeriesOf(plot);
        Assert.That(series.Count, Is.GreaterThan(0));

        var byCategory = series
            .SelectMany(s => s.Items.OfType<StackedColumnItem>().Select(i => (Series: s, Item: i)))
            .GroupBy(x => x.Item.CategoryIndex);

        foreach (var category in byCategory)
        {
            // Preserve the order the series were added in, which is the stacking order.
            var ordered = category.OrderBy(x => series.IndexOf(x.Series)).ToList();
            double expectedFloor = ordered[0].Series.BaseValue;

            foreach (var (_, item) in ordered)
            {
                Assert.That(item.Y0, Is.EqualTo(expectedFloor).Within(1e-9),
                    $"category {item.CategoryIndex}: segment does not start on the one below it");
                Assert.That(item.Y1, Is.GreaterThanOrEqualTo(item.Y0),
                    $"category {item.CategoryIndex}: segment has negative height");
                expectedFloor = item.Y1;
            }
        }
    }

    /// <summary>
    /// Under a log y-axis the stack has to start at 0.1 rather than 0, since a log axis cannot
    /// render 0. ColumnSeries.BaseValue did this; StackedColumnSeries has to place it itself.
    /// </summary>
    [Test]
    public void LogScaleStacksFromPointOneRatherThanZero()
    {
        var parameters = new PlotModelStatParameters
        {
            GroupingProperty = "None",
            MinRelativeCutoff = 0,
            MaxRelativeCutoff = 100,
            AllowAmbiguousGroups = true,
            NormalizeHistogramToFile = false,
            UseLogScaleYAxis = true
        };

        var plot = BuildPlot("Histogram of Precursor Charges", parameters);
        var series = SeriesOf(plot);
        Assert.That(series, Is.Not.Empty);
        Assert.That(series.All(s => s.BaseValue == 0.1), Is.True);

        var first = series.First();
        Assert.That(first.Items.OfType<StackedColumnItem>().All(i => i.Y0 >= 0.1 - 1e-9), Is.True,
            "a bar starts below the log axis minimum");
    }

    /// <summary>
    /// The numbers themselves: total bar height per category must equal the bin count
    /// SpectrumMatchHistograms produced. This is the check that survives any future renderer swap,
    /// and it is why the binning was extracted in phase 1.
    /// </summary>
    [TestCase("Histogram of Precursor Charges", 3)]
    [TestCase("Histogram of Fragment Charges", 4)]
    [TestCase("Histogram of Missed Cleavages", 9)]
    [TestCase("Histogram of Notch (Ambiguous PSMs Split Across Notches)", 14)]
    public void StackedBarHeightsEqualTheExtractedBins(string plotName, int plotType)
    {
        var parameters = new PlotModelStatParameters
        {
            GroupingProperty = "None",
            MinRelativeCutoff = 0,
            MaxRelativeCutoff = 100,
            AllowAmbiguousGroups = true,
            NormalizeHistogramToFile = false,
            UseLogScaleYAxis = false
        };

        MetaDrawSettings.ResetSettings();
        var bySourceFile = LoadFixture();
        var all = new ObservableCollection<SpectrumMatchFromTsv>(bySourceFile.Values.SelectMany(v => v));
        var plot = new PlotModelStat(plotName, all, bySourceFile, parameters);

        var expected = SpectrumMatchHistograms.GetHistogramData(plotType, bySourceFile);
        double expectedTotal = expected.DictsBySourceFile.Values.Sum(d => d.Values.Sum());

        double renderedTotal = SeriesOf(plot)
            .SelectMany(s => s.Items.OfType<StackedColumnItem>())
            .Sum(i => i.Value);

        Assert.That(renderedTotal, Is.EqualTo(expectedTotal).Within(1e-6),
            $"{plotName}: the bars do not add up to the bins");

        // And the per-category split, not just the total.
        var perCategory = SeriesOf(plot)
            .SelectMany(s => s.Items.OfType<StackedColumnItem>())
            .GroupBy(i => i.bin)
            .ToDictionary(g => g.Key, g => g.Sum(i => i.Value));

        foreach (var dict in expected.DictsBySourceFile.Values)
        {
            foreach (var (bin, count) in dict)
            {
                if (!perCategory.ContainsKey(bin))
                    continue; // relative-cutoff filtering can drop a bin; totals above still hold
                Assert.That(perCategory[bin], Is.GreaterThanOrEqualTo(count - 1e-6),
                    $"{plotName}: bin {bin} lost counts between binning and rendering");
            }
        }
    }

    /// <summary>
    /// Hovering a bar still names its bin, value and total. The tracker strings used {2} for the
    /// value, which was BarSeriesBase's third positional argument; RectangleBarSeries passes
    /// different arguments, so the value is now bound by name and this is the check that it resolves.
    /// </summary>
    [Test]
    public void TrackerTextNamesTheBinValueAndTotal()
    {
        var plot = BuildPlot("Histogram of Precursor Charges");
        Render(plot.Model);

        var xAxis = plot.Model.Axes.First(a => a.Position == AxisPosition.Bottom);
        var yAxis = plot.Model.Axes.First(a => a.Position == AxisPosition.Left);
        var series = SeriesOf(plot).First();
        var item = series.Items.OfType<StackedColumnItem>().First(i => i.Y1 > i.Y0);

        // Aim at the middle of the bar.
        var point = new ScreenPoint(
            (xAxis.Transform(item.X0) + xAxis.Transform(item.X1)) / 2,
            (yAxis.Transform(item.Y0) + yAxis.Transform(item.Y1)) / 2);

        var hit = series.GetNearestPoint(point, false);
        Assert.That(hit, Is.Not.Null, "hovering the middle of a bar hit nothing");
        Assert.That(hit.Text, Does.Contain($"Bin: {item.bin}"));
        Assert.That(hit.Text, Does.Contain($"Total: {item.total}"));
        Assert.That(hit.Text, Does.Contain($"{series.Title}: {item.Value}"),
            $"the tracker does not report the bar's value. Text was:\n{hit.Text}");
        Assert.That(hit.Text, Does.Not.Contain("{"), $"an unresolved placeholder remains: {hit.Text}");
    }

    private static double Round(double value) => System.Math.Round(value, 6);
}
