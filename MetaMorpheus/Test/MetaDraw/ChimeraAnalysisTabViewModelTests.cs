using NUnit.Framework;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Windows.Controls;
using System.Windows.Input;
using GuiFunctions;
using EngineLayer;
using MassSpectrometry;
using OxyPlot.Wpf;

namespace Test.MetaDraw;

[TestFixture, Apartment(System.Threading.ApartmentState.STA)]
public class ChimeraAnalysisTabViewModelTests
{
    [SetUp]
    public void SetUp()
    {
        MessageBoxHelper.SuppressMessageBoxes = true;
    }

    [Test]
    public void Constructor_InitializesProperties()
    {
        // Arrange
        var allPsms = ChimeraPlottingTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraPlottingTests.DataFile } };
        string exportDir = Path.Combine(Path.GetTempPath(), "ChimeraAnalysisTabViewModelTests");

        // Act
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles, exportDir);

        // Assert
        Assert.That(vm.ChimeraLegendViewModel, Is.Not.Null);
        Assert.That(vm.ChimeraGroupViewModels, Is.Not.Null.And.Not.Empty);
        Assert.That(vm.ExportDirectory, Is.EqualTo(exportDir));
        Assert.That(vm.SelectedExportType, Is.EqualTo("Png"));
        Assert.That(vm.ExportTypes, Is.EquivalentTo(new[] { "Pdf", "Png", "Svg" }));
        Assert.That(vm.ExportMs1Command, Is.Not.Null);
        Assert.That(vm.ExportMs2Command, Is.Not.Null);
        Assert.That(vm.ExportSequenceCoverageCommand, Is.Not.Null);
        Assert.That(vm.ExportLegendCommand, Is.Not.Null);
    }

    [Test]
    public void SelectedChimeraGroup_Setter_UpdatesLegendItemsAndNotifies()
    {
        // Arrange
        var allPsms = ChimeraPlottingTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraPlottingTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles);
        var group = vm.ChimeraGroupViewModels[0];
        bool propertyChanged = false;
        vm.PropertyChanged += (s, e) => { if (e.PropertyName == nameof(vm.SelectedChimeraGroup)) propertyChanged = true; };

        // Act
        vm.SelectedChimeraGroup = group;

        // Assert
        Assert.That(vm.SelectedChimeraGroup, Is.EqualTo(group));
        Assert.That(vm.ChimeraLegendViewModel.ChimeraLegendItems, Is.EqualTo(group.LegendItems));
        Assert.That(propertyChanged, Is.True);
    }

    [Test]
    public void ChimeraLegendViewModel_Setter_Notifies()
    {
        // Arrange
        var allPsms = ChimeraPlottingTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraPlottingTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles);
        var legendVm = new ChimeraLegendViewModel();
        bool propertyChanged = false;
        vm.PropertyChanged += (s, e) => { if (e.PropertyName == nameof(vm.ChimeraLegendViewModel)) propertyChanged = true; };

        // Act
        vm.ChimeraLegendViewModel = legendVm;

        // Assert
        Assert.That(vm.ChimeraLegendViewModel, Is.EqualTo(legendVm));
        Assert.That(propertyChanged, Is.True);
    }

    [Test]
    public void ChimeraDrawnSequence_Setter_Notifies()
    {
        // Arrange
        var allPsms = ChimeraPlottingTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraPlottingTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles);
        var group = vm.ChimeraGroupViewModels[0];
        var canvas = new Canvas();
        var drawnSeq = new ChimeraDrawnSequence(canvas, group);
        bool propertyChanged = false;
        vm.PropertyChanged += (s, e) => { if (e.PropertyName == nameof(vm.ChimeraDrawnSequence)) propertyChanged = true; };

        // Act
        vm.ChimeraDrawnSequence = drawnSeq;

        // Assert
        Assert.That(vm.ChimeraDrawnSequence, Is.EqualTo(drawnSeq));
        Assert.That(propertyChanged, Is.True);
    }

    [Test]
    public void UseLetterOnly_Setter_UpdatesIonColorsAndNotifies()
    {
        // Arrange
        var allPsms = ChimeraPlottingTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraPlottingTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles);
        bool propertyChanged = false;
        vm.PropertyChanged += (s, e) => { if (e.PropertyName == nameof(vm.UseLetterOnly)) propertyChanged = true; };

        // Act
        vm.UseLetterOnly = true;

        // Assert
        Assert.That(vm.UseLetterOnly, Is.True);
        Assert.That(propertyChanged, Is.True);
        // Check that AssignIonColors was called by verifying IsColorInitialized is true for all groups
        Assert.That(vm.ChimeraGroupViewModels.All(g => g.MatchedFragmentIonsByColor != null), Is.True);
    }

    [Test]
    public void ExportDirectory_Setter_CreatesDirectoryAndNotifies()
    {
        // Arrange
        var allPsms = ChimeraPlottingTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraPlottingTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles);
        string tempDir = Path.Combine(Path.GetTempPath(), "ChimeraAnalysisTabViewModelTests_ExportDir");
        if (Directory.Exists(tempDir))
            Directory.Delete(tempDir, true);
        bool propertyChanged = false;
        vm.PropertyChanged += (s, e) => { if (e.PropertyName == nameof(vm.ExportDirectory)) propertyChanged = true; };

        // Act
        vm.ExportDirectory = tempDir;

        // Assert
        Assert.That(vm.ExportDirectory, Is.EqualTo(tempDir));
        Assert.That(Directory.Exists(tempDir), Is.True);
        Assert.That(propertyChanged, Is.True);

        // Cleanup
        Directory.Delete(tempDir, true);
    }

    [Test]
    public void SelectedExportType_Setter_Notifies()
    {
        // Arrange
        var allPsms = ChimeraPlottingTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraPlottingTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles);
        bool propertyChanged = false;
        vm.PropertyChanged += (s, e) => { if (e.PropertyName == nameof(vm.SelectedExportType)) propertyChanged = true; };

        // Act
        vm.SelectedExportType = "Pdf";

        // Assert
        Assert.That(vm.SelectedExportType, Is.EqualTo("Pdf"));
        Assert.That(propertyChanged, Is.True);
    }

    [Test]
    public void ExportMs1Command_DoesNothing_WhenSelectedChimeraGroupIsNull()
    {
        // Arrange
        var allPsms = ChimeraPlottingTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraPlottingTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles);

        string tempDir = Path.Combine(Path.GetTempPath(), "ChimeraAnalysisTabViewModelTests_ExportMs1_Null");
        vm.ExportDirectory = tempDir;
        Directory.CreateDirectory(tempDir);

        // Act & Assert
        Assert.DoesNotThrow(() => vm.ExportMs1Command.Execute(null));
        Assert.That(Directory.GetFiles(tempDir, "*_MS1.*").Length, Is.EqualTo(0));

        // Cleanup
        Directory.Delete(tempDir, true);
    }

    [Test]
    public void ExportMs2Command_DoesNothing_WhenSelectedChimeraGroupIsNull()
    {
        // Arrange
        var allPsms = ChimeraPlottingTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraPlottingTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles);

        string tempDir = Path.Combine(Path.GetTempPath(), "ChimeraAnalysisTabViewModelTests_ExportMs2_Null");
        vm.ExportDirectory = tempDir;
        Directory.CreateDirectory(tempDir);

        // Act & Assert
        Assert.DoesNotThrow(() => vm.ExportMs2Command.Execute(null));
        Assert.That(Directory.GetFiles(tempDir, "*_MS2.*").Length, Is.EqualTo(0));

        // Cleanup
        Directory.Delete(tempDir, true);
    }

    [Test]
    public void ExportLegendCommand_DoesNothing_WhenArgumentIsNull()
    {
        // Arrange
        var allPsms = ChimeraPlottingTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraPlottingTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles);
        string tempDir = Path.Combine(Path.GetTempPath(), "ChimeraAnalysisTabViewModelTests_ExportLegend_Null");
        vm.ExportDirectory = tempDir;
        Directory.CreateDirectory(tempDir);

        // Act & Assert
        Assert.DoesNotThrow(() => vm.ExportLegendCommand.Execute(null));
        Assert.That(Directory.GetFiles(tempDir, "*_Legend.*").Length, Is.EqualTo(0));

        // Cleanup
        Directory.Delete(tempDir, true);
    }

    [TestCase("Pdf")]
    [TestCase("Png")]
    [TestCase("Svg")]
    public void ExportMs1Command_ExportsCorrectFormat(string format)
    {
        // Arrange
        var allPsms = ChimeraPlottingTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraPlottingTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles);
        vm.SelectedChimeraGroup = vm.ChimeraGroupViewModels[0];
        vm.SelectedExportType = format;
        string tempDir = Path.Combine(Path.GetTempPath(), $"ChimeraAnalysisTabViewModelTests_ExportMs1_{format}");
        vm.ExportDirectory = tempDir;
        Directory.CreateDirectory(tempDir);
        vm.Ms1ChimeraPlot = new Ms1ChimeraPlot(new OxyPlot.Wpf.PlotView(), vm.SelectedChimeraGroup);

        // Act
        Assert.DoesNotThrow(() => vm.ExportMs1Command.Execute(null));

        // Assert
        string ext = format.ToLower();
        Assert.That(Directory.GetFiles(tempDir, $"*_MS1.{ext}").Length, Is.GreaterThan(0));

        // Cleanup
        foreach (var file in Directory.GetFiles(tempDir))
            File.Delete(file);
        Directory.Delete(tempDir, true);
    }

    [TestCase("Pdf")]
    [TestCase("Png")]
    [TestCase("Svg")]
    public void ExportMs2Command_ExportsCorrectFormat(string format)
    {
        // Arrange
        var allPsms = ChimeraPlottingTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraPlottingTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles);
        vm.SelectedChimeraGroup = vm.ChimeraGroupViewModels[0];
        vm.SelectedExportType = format;
        string tempDir = Path.Combine(Path.GetTempPath(), $"ChimeraAnalysisTabViewModelTests_ExportMs2_{format}");
        vm.ExportDirectory = tempDir;
        Directory.CreateDirectory(tempDir);
        var plotView = new OxyPlot.Wpf.PlotView { Width = 200, Height = 100 };
        vm.ChimeraSpectrumMatchPlot = new ChimeraSpectrumMatchPlot(plotView, vm.SelectedChimeraGroup);

        // Act
        Assert.DoesNotThrow(() => vm.ExportMs2Command.Execute(null));

        // Assert
        string ext = format.ToLower();
        Assert.That(Directory.GetFiles(tempDir, $"*_MS2.{ext}").Length, Is.GreaterThan(0));

        // Cleanup
        foreach (var file in Directory.GetFiles(tempDir))
            File.Delete(file);
        Directory.Delete(tempDir, true);
    }

    [Test]
    public void ExportSequenceCoverageCommand_InvokesExportSequenceCoverage()
    {
        // Arrange
        var allPsms = ChimeraPlottingTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraPlottingTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles);
        vm.SelectedChimeraGroup = vm.ChimeraGroupViewModels[0];
        vm.SelectedExportType = "Png";
        string tempDir = Path.Combine(Path.GetTempPath(), "ChimeraAnalysisTabViewModelTests_ExportSeqCov");
        vm.ExportDirectory = tempDir;
        Directory.CreateDirectory(tempDir);

        // Setup a drawn sequence for export
        var canvas = new Canvas { Width = 200, Height = 100 };
        // Set actual width and height to match the intended export size
        canvas.Children.Add(new System.Windows.Shapes.Rectangle
        {
            Width = 200,
            Height = 100,
            Fill = System.Windows.Media.Brushes.Transparent
        });
        canvas.Measure(new System.Windows.Size(200, 100));
        canvas.Arrange(new System.Windows.Rect(0, 0, 200, 100));
        canvas.UpdateLayout();

        vm.ChimeraDrawnSequence = new ChimeraDrawnSequence(canvas, vm.SelectedChimeraGroup);

        // Act
        Assert.DoesNotThrow(() => vm.ExportSequenceCoverageCommand.Execute(new()));

        // Assert
        string[] files = Directory.GetFiles(tempDir, "*_SequenceCoverage.png");
        Assert.That(files.Length, Is.GreaterThan(0));

        // Cleanup
        foreach (var file in files)
            File.Delete(file);
        Directory.Delete(tempDir, true);
    }

    [Test]
    public void ExportLegendCommand_InvokesExportLegend()
    {
        // Arrange
        var allPsms = ChimeraPlottingTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraPlottingTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles);
        vm.SelectedChimeraGroup = vm.ChimeraGroupViewModels[0];
        vm.SelectedExportType = "Png";
        string tempDir = Path.Combine(Path.GetTempPath(), "ChimeraAnalysisTabViewModelTests_ExportLegend");
        vm.ExportDirectory = tempDir;
        Directory.CreateDirectory(tempDir);

        // Setup a dummy legend FrameworkElement
        var legend = new Canvas { Width = 100, Height = 50 };
        // Add a child so GetDescendantBounds returns a valid rect
        legend.Children.Add(new System.Windows.Shapes.Rectangle
        {
            Width = 100,
            Height = 50,
            Fill = System.Windows.Media.Brushes.Transparent
        });
        legend.Measure(new System.Windows.Size(100, 50));
        legend.Arrange(new System.Windows.Rect(0, 0, 100, 50));
        legend.UpdateLayout();
        legend.SetValue(Canvas.WidthProperty, 100.0);
        legend.SetValue(Canvas.HeightProperty, 50.0);

        // Act
        Assert.DoesNotThrow(() => vm.ExportLegendCommand.Execute(legend));

        // Assert
        string[] files = Directory.GetFiles(tempDir, "*_Legend.png");
        Assert.That(files.Length, Is.GreaterThan(0));

        // Cleanup
        foreach (var file in files)
            File.Delete(file);
        Directory.Delete(tempDir, true);
    }
}
