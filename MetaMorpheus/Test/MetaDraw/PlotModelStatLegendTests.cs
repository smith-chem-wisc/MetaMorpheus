using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using GuiFunctions;
using NUnit.Framework;
using OxyPlot;
using Readers;

namespace Test.MetaDraw;

/// <summary>
/// OxyPlot 2.0 drew a legend on every PlotModel; 2.1 moved legends into PlotModel.Legends and draws
/// none when that is empty. The scatter plots add title-only series purely as legend keys, so losing
/// the legend there is silent. These render to SVG and look for the legend text.
/// </summary>
public class PlotModelStatLegendTests
{
    private static string RenderSvg(string plotName)
    {
        MetaDrawSettings.ResetSettings();
        string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "VariantCrossTest.psmtsv");
        var matches = SpectrumMatchTsvReader.ReadTsv(path, out _);
        var bySourceFile = matches.GroupBy(p => p.FileNameWithoutExtension)
            .ToDictionary(g => g.Key, g => new ObservableCollection<SpectrumMatchFromTsv>(g));
        var plot = new PlotModelStat(plotName, new ObservableCollection<SpectrumMatchFromTsv>(matches), bySourceFile);
        return new SvgExporter { Width = 800, Height = 600, IsDocument = true }.ExportToString(plot.Model);
    }

    [TestCase("Precursor PPM Error vs. RT")]
    [TestCase("Predicted RT vs. Observed RT")]
    public void ScatterPlotsDrawTheVariantKey(string plotName)
    {
        // VariantCrossTest.psmtsv has both kinds, so both keys are added.
        string svg = RenderSvg(plotName);
        Assert.That(svg, Does.Contain(">non-variant PSMs<"));
        Assert.That(svg, Does.Contain(">variant PSMs<"));
    }

    [Test]
    public void HistogramsDrawTheSourceFileLegendTitle()
    {
        Assert.That(RenderSvg("Histogram of Precursor Charges"), Does.Contain(">Source file(s)<"));
    }
}
