using OxyPlot;
using OxyPlot.Wpf;
using System.Linq;
using NUnit.Framework;
using GuiFunctions;
using System.Windows.Controls;
using System.IO;
using System.Threading;

namespace Test.MetaDraw;

[TestFixture, Apartment(ApartmentState.STA)]
public class ChimeraSpectrumMatchPlotTests
{
    [SetUp]
    public void SetUp()
    {
        MessageBoxHelper.SuppressMessageBoxes = true;
    }

    [Test]
    public void ChimeraSpectrumMatchPlot_Constructor_InitializesProperties()
    {
        // Arrange
        var chimeraGroup = ChimeraPlottingTests.OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        var plotView = new PlotView();

        // Act
        var plot = new ChimeraSpectrumMatchPlot(plotView, chimeraGroup);

        // Assert
        Assert.That(plot.Model, Is.Not.Null);
        Assert.That(plot.MatchedFragmentIons, Is.Not.Null.And.Not.Empty, "MatchedFragmentIons should be initialized and not empty.");
        Assert.That(plot.Model.Series, Is.Not.Empty, "Model should have at least one series after construction.");
    }

    [Test]
    public void ChimeraSpectrumMatchPlot_AnnotatesMatchedIonsWithCorrectColors()
    {
        // Arrange
        var chimeraGroup = ChimeraPlottingTests.TwoProteinsTwoProteoformChimeraGroup.ChimeraGroup;
        var expectedColors = chimeraGroup.MatchedFragmentIonsByColor.Keys;
        var plotView = new PlotView();

        // Act
        var plot = new ChimeraSpectrumMatchPlot(plotView, chimeraGroup);

        // Assert
        // For each expected color, there should be at least one series with that color
        foreach (var color in expectedColors)
        {
            bool found = plot.Model.Series
                .OfType<OxyPlot.Series.LineSeries>()
                .Any(ls => ls.Color == color);
            Assert.That(found, Is.True, $"Expected to find a peak with color {color} for matched fragment ion.");
        }
    }

    [Test]
    public void ChimeraSpectrumMatchPlot_SetsMzMaxIfProvided()
    {
        // Arrange
        var chimeraGroup = ChimeraPlottingTests.OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        var plotView = new PlotView();
        double mzMax = 1000.0;

        // Act
        var plot = new ChimeraSpectrumMatchPlot(plotView, chimeraGroup, mzMax);

        // Assert
        Assert.That(plot.Model.Axes[0].Maximum, Is.EqualTo(mzMax).Within(1e-6));
    }

    [Test]
    public void ChimeraSpectrumMatchPlot_ExportPlot_CreatesExportedFile()
    {
        // Arrange
        var chimeraGroup = ChimeraPlottingTests.OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        var plotView = new PlotView();
        var plot = new ChimeraSpectrumMatchPlot(plotView, chimeraGroup);
        string tempDir = Path.Combine(Path.GetTempPath(), "ChimeraSpectrumMatchPlotTests");
        Directory.CreateDirectory(tempDir);
        string exportPath = Path.Combine(tempDir, "test_export.png");

        // Act
        plot.ExportPlot(exportPath);

        // Assert
        Assert.That(File.Exists(exportPath), Is.True, "Exported plot file should exist.");

        // Cleanup
        File.Delete(exportPath);
        Directory.Delete(tempDir);
    }

    [Test]
    public void ChimeraSpectrumMatchPlot_ExportPlot_WithLegend_CreatesExportedFile()
    {
        // Arrange
        var chimeraGroup = ChimeraPlottingTests.OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        var plotView = new PlotView();
        var plot = new ChimeraSpectrumMatchPlot(plotView, chimeraGroup);

        // Create a simple legend canvas
        var legend = new Canvas
        {
            Width = 100,
            Height = 50
        };
        // Set actual width and height to match the intended export size
        legend.Measure(new System.Windows.Size(100, 50));
        legend.Arrange(new System.Windows.Rect(0, 0, 100, 50));
        legend.UpdateLayout();
        legend.SetValue(Canvas.WidthProperty, 100.0);
        legend.SetValue(Canvas.HeightProperty, 50.0);

        // Ensure legend is shown
        bool originalShowLegend = MetaDrawSettings.ShowLegend;
        MetaDrawSettings.ShowLegend = true;

        string tempDir = Path.Combine(Path.GetTempPath(), "ChimeraSpectrumMatchPlotTests_Legend");
        Directory.CreateDirectory(tempDir);
        string exportPath = Path.Combine(tempDir, "test_export_with_legend.png");

        try
        {
            // Act
            plot.ExportPlot(exportPath, legend);

            // Assert
            Assert.That(File.Exists(exportPath), Is.True, "Exported plot file with legend should exist.");
        }
        finally
        {
            // Cleanup
            if (File.Exists(exportPath))
                File.Delete(exportPath);
            if (Directory.Exists(tempDir))
                Directory.Delete(tempDir);
            MetaDrawSettings.ShowLegend = originalShowLegend;
        }
    }


    [Test]
    public void ChimeraSpectrumMatchPlot_StaticColorDictionaries_AreInitialized()
    {
        // Act & Assert
        Assert.That(ChimeraSpectrumMatchPlot.MultipleProteinSharedColor, Is.EqualTo(OxyColors.Black));
        Assert.That(ChimeraSpectrumMatchPlot.ColorByProteinDictionary, Is.Not.Null.And.Not.Empty);
        Assert.That(ChimeraSpectrumMatchPlot.OverflowColors, Is.Not.Null.And.Not.Empty);
    }
}
