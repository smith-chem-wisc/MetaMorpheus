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
using System.Drawing.Imaging;
using System;
using System.Drawing;
using System.Reflection;
using GuiFunctions.MetaDraw;
using System.Windows;
using Point = System.Windows.Point;
using Size = System.Windows.Size;
using Readers;

namespace Test.MetaDraw;

[TestFixture, Apartment(System.Threading.ApartmentState.STA)]
public class ChimeraAnalysisTabViewModelTests
{
    public static string TestExportDirectory => Path.Combine(TestContext.CurrentContext.TestDirectory, "MetaDraw", "ChimeraAnalysisTabViewModelTests");

    [OneTimeSetUp]
    public static void OneTimeSetup()
    {
        MessageBoxHelper.SuppressMessageBoxes = true;
        GlobalVariables.AnalyteType = AnalyteType.Proteoform;
        // Ensure the export directory exists in a new state
        if (Directory.Exists(TestExportDirectory))
        {
            Directory.Delete(TestExportDirectory, true);
        }
        Directory.CreateDirectory(TestExportDirectory);
    }

    [OneTimeTearDown]
    public static void OneTimeTearDown()
    {
        GlobalVariables.AnalyteType = AnalyteType.Peptide;
        // Clean up the export directory after tests
        if (Directory.Exists(TestExportDirectory))
        {
            Directory.Delete(TestExportDirectory, true);
        }
    }

    [Test]
    [NonParallelizable]
    public void ChimeraTabViewModel_ResultsFilteredCorrectly()
    {
        var dataFiles = new Dictionary<string, MsDataFile>()
        {
            {"FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile }
        };

        var chimeraAnalysisTabViewModel = new ChimeraAnalysisTabViewModel(TestExportDirectory);
        chimeraAnalysisTabViewModel.ProcessChimeraData(ChimeraGroupViewModelTests.AllMatches, dataFiles);
        Assert.That(chimeraAnalysisTabViewModel.ChimeraGroupViewModels.All(p => p.Count > 1),
            Is.True, "All chimera groups should have at least two PSMs.");
        Assert.That(chimeraAnalysisTabViewModel.ChimeraGroupViewModels.Count, Is.GreaterThan(0), "There should be at least one chimera group.");

        // Check that all PSMs in chimera groups pass the Q-value filter and do not contain decoys if ShowDecoys is false
        chimeraAnalysisTabViewModel.ChimeraGroupViewModels.SelectMany(p => p.ChimericPsms)
            .ToList()
            .ForEach(p =>
            {
                Assert.That(p.Psm.QValue, Is.LessThanOrEqualTo(MetaDrawSettings.QValueFilter), "Chimeric PSMs should pass the Q-value filter.");
                Assert.That(p.Psm.DecoyContamTarget, Does.Not.Contain('D'), "Chimeric PSMs should not contain decoys if ShowDecoys is false.");
            });

        MetaDrawSettings.QValueFilter = 0.5; // set filter higher
        MetaDrawSettings.ShowDecoys = true; // show decoys

        chimeraAnalysisTabViewModel = new ChimeraAnalysisTabViewModel(TestExportDirectory);
        chimeraAnalysisTabViewModel.ProcessChimeraData(ChimeraGroupViewModelTests.AllMatches, dataFiles);
        Assert.That(chimeraAnalysisTabViewModel.ChimeraGroupViewModels.All(p => p.Count > 1),
            Is.True, "All chimera groups should have at least two PSMs after changing settings.");
        Assert.That(chimeraAnalysisTabViewModel.ChimeraGroupViewModels.Count, Is.GreaterThan(0), "There should still be at least one chimera group after changing settings.");

        // Check that all PSMs in chimera groups pass the new Q-value filter and may contain decoys
        chimeraAnalysisTabViewModel.ChimeraGroupViewModels.SelectMany(p => p.ChimericPsms)
            .ToList()
            .ForEach(p =>
            {
                Assert.That(p.Psm.QValue, Is.LessThanOrEqualTo(MetaDrawSettings.QValueFilter), "Chimeric PSMs should pass the new Q-value filter.");
            });
        Assert.That(chimeraAnalysisTabViewModel.ChimeraGroupViewModels.Any(p => p.ChimericPsms.Any(q => q.Psm.DecoyContamTarget.Contains('D'))),
            Is.True, "Some chimeric PSMs should contain decoys after changing ShowDecoys to true.");

        // Reset settings for other tests
        MetaDrawSettings.QValueFilter = 0.01; // reset filter
        MetaDrawSettings.ShowDecoys = false; // reset to not show decoys
    }

    [Test]
    [NonParallelizable]
    public static void ChimeraTabViewModel_PrecursorAssignmentIsCorrectInGroups()
    {
        var dataFiles = new Dictionary<string, MsDataFile>()
        {
            {"FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile }
        };
        var chimeraAnalysisTabViewModel = new ChimeraAnalysisTabViewModel(TestExportDirectory);
        chimeraAnalysisTabViewModel.ProcessChimeraData(ChimeraGroupViewModelTests.AllMatches, dataFiles);

        foreach (var chimeraGroup in chimeraAnalysisTabViewModel.ChimeraGroupViewModels)
        {
            // Only test groups with at least 2 PSMs
            if (chimeraGroup.Count < 2)
                continue;

            var chimericPsms = chimeraGroup.ChimericPsms.ToList();
            for (int i = 0; i < chimericPsms.Count; i++)
            {
                var psm = chimericPsms[i].Psm;
                var envelope = chimericPsms[i].PrecursorEnvelope;

                // Compare to all other envelopes in the group
                double minDiff = Math.Abs(psm.PrecursorMass - envelope.MonoisotopicMass);
                for (int j = 0; j < chimericPsms.Count; j++)
                {
                    if (i == j) continue;

                    var otherEnvelope = chimericPsms[j].PrecursorEnvelope;
                    if (psm.ChargeState != otherEnvelope.Charge) 
                        continue;

                    double diff = Math.Abs(psm.PrecursorMass - otherEnvelope.MonoisotopicMass);
                    Assert.That(minDiff, Is.LessThanOrEqualTo(diff),
                        $"PSM {i} in group (scan {chimeraGroup.Ms2ScanNumber}) should be assigned to its closest envelope.");
                }
            }
        }
    }

    [Test]
    public void Constructor_InitializesProperties()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        string exportDir = Path.Combine(Path.GetTempPath(), "ChimeraAnalysisTabViewModelTests");

        // Act
        var vm = new ChimeraAnalysisTabViewModel(exportDir);
        vm.ProcessChimeraData(allPsms, dataFiles);

        // Assert
        Assert.That(vm.ChimeraGroupViewModels, Is.Not.Null.And.Not.Empty);
        Assert.That(vm.ExportDirectory, Is.EqualTo(exportDir));
        Assert.That(vm.ExportMs1Command, Is.Not.Null);
        Assert.That(vm.ExportMs2Command, Is.Not.Null);
        Assert.That(vm.ExportSequenceCoverageCommand, Is.Not.Null);
        Assert.That(vm.ExportLegendCommand, Is.Not.Null);
    }

    [Test]
    public void ChimeraDrawnSequence_Setter_Notifies()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);
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
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);
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
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);
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
    public void ExportMs1Command_DoesNothing_WhenSelectedChimeraGroupIsNull()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);

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
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);

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
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);
        string tempDir = Path.Combine(Path.GetTempPath(), "ChimeraAnalysisTabViewModelTests_ExportLegend_Null");
        vm.ExportDirectory = tempDir;
        Directory.CreateDirectory(tempDir);

        // Act & Assert
        Assert.DoesNotThrow(() => vm.ExportLegendCommand.Execute(null));
        Assert.That(Directory.GetFiles(tempDir, "*_Legend.*").Length, Is.EqualTo(0));

        // Cleanup
        Directory.Delete(tempDir, true);
    }

    [Test]
    public void ExportMs1Command_ExportsCorrectFormat()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);
        vm.SelectedChimeraGroup = vm.ChimeraGroupViewModels[0];
        string tempDir = Path.Combine(Path.GetTempPath(), $"ChimeraAnalysisTabViewModelTests_ExportMs1_{MetaDrawSettings.ExportType}");
        vm.ExportDirectory = tempDir;
        Directory.CreateDirectory(tempDir);

        var plotView = new OxyPlot.Wpf.PlotView { Width = 200, Height = 100 };
        plotView.Measure(new Size(200, 100));
        plotView.Arrange(new Rect(0, 0, 200, 100));
        plotView.UpdateLayout();

        vm.Ms1ChimeraPlot = new Ms1ChimeraPlot(plotView, vm.SelectedChimeraGroup);

        // Act
        Assert.DoesNotThrow(() => vm.ExportMs1Command.Execute(null));

        // Assert
        string ext = MetaDrawSettings.ExportType.ToLower();
        Assert.That(Directory.GetFiles(tempDir, $"*_MS1.{ext}").Length, Is.GreaterThan(0));

        // Cleanup
        foreach (var file in Directory.GetFiles(tempDir))
            File.Delete(file);
        Directory.Delete(tempDir, true);
    }

    [Test]
    public void ExportMs2Command_ExportsCorrectFormat()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);
        vm.SelectedChimeraGroup = vm.ChimeraGroupViewModels[0];
        string tempDir = Path.Combine(Path.GetTempPath(), $"ChimeraAnalysisTabViewModelTests_ExportMs2_{MetaDrawSettings.ExportType}");
        vm.ExportDirectory = tempDir;
        Directory.CreateDirectory(tempDir);

        var plotView = new OxyPlot.Wpf.PlotView { Width = 200, Height = 100 };
        plotView.Measure(new Size(200, 100));
        plotView.Arrange(new Rect(0, 0, 200, 100));
        plotView.UpdateLayout();

        vm.ChimeraSpectrumMatchPlot = new ChimeraSpectrumMatchPlot(plotView, vm.SelectedChimeraGroup);

        // Act
        Assert.DoesNotThrow(() => vm.ExportMs2Command.Execute(null));

        // Assert
        string ext = MetaDrawSettings.ExportType.ToLower();
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
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);
        vm.SelectedChimeraGroup = vm.ChimeraGroupViewModels[0];
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
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);
        vm.SelectedChimeraGroup = vm.ChimeraGroupViewModels[0];
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
        string[] files = Directory.GetFiles(tempDir, "*_Legend.pdf");
        Assert.That(files.Length, Is.GreaterThan(0));

        // Cleanup
        foreach (var file in files)
            File.Delete(file);
        Directory.Delete(tempDir, true);
    }

    [Test]
    public void ConstructChimericPsms_ContinuesIfDataFileNotFound()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        // Provide an empty dataFiles dictionary so no data file can be found
        var dataFiles = new Dictionary<string, MsDataFile>();

        // Act
        var result = ChimeraAnalysisTabViewModel.ConstructChimericPsms(allPsms, dataFiles);

        // Assert
        Assert.That(result, Is.Empty, "No groups should be constructed if no data file is found.");
    }

    [Test]
    public void PrecursorIonAnnotations_AreLetterOnly_WhenUseLetterOnlyIsTrue()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);

        // Act
        vm.UseLetterOnly = true;

        // Assert
        foreach (var group in vm.ChimeraGroupViewModels)
        {
            foreach (var kvp in group.PrecursorIonsByColor)
            {
                foreach (var tuple in kvp.Value)
                {
                    if (!tuple.Item2.Contains("Miso"))
                    {
                        Assert.That(tuple.Item2, Is.Not.Empty, "Annotation should not be empty when using letter only.");
                        Assert.That(tuple.Item2.Length, Is.EqualTo(1), "Annotation should be a single letter when using letter only.");
                        Assert.That("ABCDEFGHIJKLMNOPQRSTUVWXYZ".Contains(tuple.Item2), Is.True, "Annotation should be a letter.");
                    }
                }
            }
        }
    }

    [Test]
    public void ChimeraSpectrumMatchPlot_GetterSetter_Works()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);
        var plotView = new OxyPlot.Wpf.PlotView();
        var group = vm.ChimeraGroupViewModels[0];
        var plot = new ChimeraSpectrumMatchPlot(plotView, group);

        // Act
        vm.ChimeraSpectrumMatchPlot = plot;

        // Assert
        Assert.That(vm.ChimeraSpectrumMatchPlot, Is.EqualTo(plot));
    }

    [Test]
    public void Ms1ChimeraPlot_GetterSetter_Works()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);
        var plotView = new OxyPlot.Wpf.PlotView();
        var group = vm.ChimeraGroupViewModels[0];
        var plot = new Ms1ChimeraPlot(plotView, group);

        // Act
        vm.Ms1ChimeraPlot = plot;

        // Assert
        Assert.That(vm.Ms1ChimeraPlot, Is.EqualTo(plot));
    }

    [Test]
    public void ExportMs2Command_ExportsWithLegend_AndLegendWithinBounds()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);
        vm.SelectedChimeraGroup = vm.ChimeraGroupViewModels[0];
        string tempDir = Path.Combine(Path.GetTempPath(), "ChimeraAnalysisTabViewModelTests_ExportMs2_Legend");
        vm.ExportDirectory = tempDir;
        Directory.CreateDirectory(tempDir);

        // Setup plot and legend
        var plotView = new OxyPlot.Wpf.PlotView { Width = 400, Height = 200 };
        plotView.Measure(new System.Windows.Size(400, 200));
        plotView.Arrange(new System.Windows.Rect(0, 0, 400, 200));
        plotView.UpdateLayout();
        vm.ChimeraSpectrumMatchPlot = new ChimeraSpectrumMatchPlot(plotView, vm.SelectedChimeraGroup);

        var legend = new ChimeraLegendCanvas(vm.SelectedChimeraGroup)
        {
            Width = 100,
            Height = 50
        };
        Canvas.SetLeft(legend, 10);
        Canvas.SetTop(legend, 10);

        var parentCanvas = new Canvas { Width = 400, Height = 200 };
        parentCanvas.Children.Add(legend);
        legend.Measure(new System.Windows.Size(100, 50));
        legend.Arrange(new System.Windows.Rect(10, 10, 100, 50));
        legend.UpdateLayout();

        vm.LegendCanvas = legend;

        Assert.That(Canvas.GetLeft(legend) + legend.Width <= plotView.Width);
        Assert.That(Canvas.GetTop(legend) + legend.Height <= plotView.Height);

        // Act
        Assert.DoesNotThrow(() => vm.ExportMs2Command.Execute(null));

        // Assert
        string ext = MetaDrawSettings.ExportType.ToLower();
        Assert.That(Directory.GetFiles(tempDir, $"*_MS2.{ext}").Length, Is.GreaterThan(0));

        // Cleanup
        foreach (var file in Directory.GetFiles(tempDir))
            File.Delete(file);
        Directory.Delete(tempDir, true);
    }

    [Test]
    public void ExportAllCommand_ExportsWithLegend_AndLegendWithinBounds()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);
        vm.SelectedChimeraGroup = vm.ChimeraGroupViewModels[0];
        string tempDir = Path.Combine(Path.GetTempPath(), "ChimeraAnalysisTabViewModelTests_ExportAll_Legend");
        vm.ExportDirectory = tempDir;
        Directory.CreateDirectory(tempDir);

        // Setup plot and legend
        var ms1PlotView = new OxyPlot.Wpf.PlotView { Width = 400, Height = 200 };
        ms1PlotView.Measure(new System.Windows.Size(400, 200));
        ms1PlotView.Arrange(new System.Windows.Rect(0, 0, 400, 200));
        ms1PlotView.UpdateLayout();
        vm.Ms1ChimeraPlot = new Ms1ChimeraPlot(ms1PlotView, vm.SelectedChimeraGroup);

        var ms2PlotView = new OxyPlot.Wpf.PlotView { Width = 400, Height = 200 };
        ms2PlotView.Measure(new System.Windows.Size(400, 200));
        ms2PlotView.Arrange(new System.Windows.Rect(0, 0, 400, 200));
        ms2PlotView.UpdateLayout();
        vm.ChimeraSpectrumMatchPlot = new ChimeraSpectrumMatchPlot(ms2PlotView, vm.SelectedChimeraGroup);

        vm.ChimeraDrawnSequence = new ChimeraDrawnSequence(new Canvas { Width = 200, Height = 100 }, vm.SelectedChimeraGroup);

        var legend = new ChimeraLegendCanvas(vm.SelectedChimeraGroup)
        {
            Width = 100,
            Height = 50
        };
        Canvas.SetLeft(legend, 10);
        Canvas.SetTop(legend, 10);

        var parentCanvas = new Canvas { Width = 400, Height = 200 };
        parentCanvas.Children.Add(legend);
        legend.Measure(new System.Windows.Size(100, 50));
        legend.Arrange(new System.Windows.Rect(10, 10, 100, 50));
        legend.UpdateLayout();

        vm.LegendCanvas = legend;

        Assert.That(Canvas.GetLeft(legend) + legend.Width <= ms1PlotView.Width);
        Assert.That(Canvas.GetTop(legend) + legend.Height <= ms1PlotView.Height);

        // Act
        Assert.DoesNotThrow(() => vm.ExportAllCommand.Execute(legend));

        // Assert
        string ext = MetaDrawSettings.ExportType.ToLower();
        Assert.That(Directory.GetFiles(tempDir, $"*_ALL.{ext}").Length, Is.GreaterThan(0));

        // Cleanup
        foreach (var file in Directory.GetFiles(tempDir))
            File.Delete(file);
        Directory.Delete(tempDir, true);
    }

    
    [Test]
    public void CombinePlotAndLegend_CombinesBitmaps_WhenLegendPresent()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);

        var plotView = new Canvas { Width = 200, Height = 100, Background = System.Windows.Media.Brushes.White };
        var legend = new Canvas { Width = 50, Height = 30, Background = System.Windows.Media.Brushes.White };
        Canvas.SetLeft(legend, 10);
        Canvas.SetTop(legend, 10);

        // Add both to a parent to simulate a shared visual tree
        var parent = new Canvas { Width = 220, Height = 140 };
        parent.Children.Add(plotView);
        parent.Children.Add(legend);

        // Add to a Window and show (hidden)
        var window = new Window
        {
            Content = parent,
            Width = 220,
            Height = 140,
            WindowStyle = WindowStyle.None,
            ShowInTaskbar = false,
            ShowActivated = false,
            Visibility = Visibility.Hidden
        };
        window.Show();

        // Force layout
        parent.UpdateLayout();
        plotView.UpdateLayout();
        legend.UpdateLayout();

        // Let dispatcher process layout/render
        System.Windows.Threading.Dispatcher.CurrentDispatcher.Invoke(
            System.Windows.Threading.DispatcherPriority.Background,
            new Action(delegate { }));

        // Use reflection to access the private method
        var method = typeof(ChimeraAnalysisTabViewModel).GetMethod("CombinePlotAndLegend", BindingFlags.NonPublic | BindingFlags.Instance);
        Assert.That(method, Is.Not.Null);

        // Act
        var bitmap = method.Invoke(vm, new object[] { plotView, legend }) as Bitmap;

        // Assert
        Assert.That(bitmap, Is.Not.Null);
        Assert.That(bitmap.Width, Is.GreaterThan(0));
        Assert.That(bitmap.Height, Is.GreaterThan(0));
        bitmap.Dispose();
    }

    [Test]
    public void CombinePlotAndLegend_ReturnsPlotBitmap_WhenLegendIsNull()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);

        var plotView = new Canvas { Width = 200, Height = 100, Background = System.Windows.Media.Brushes.White };

        // Add both to a parent to simulate a shared visual tree
        var parent = new Canvas { Width = 220, Height = 140 };
        parent.Children.Add(plotView);

        // Add to a Window and show (hidden)
        var window = new Window
        {
            Content = parent,
            Width = 220,
            Height = 140,
            WindowStyle = WindowStyle.None,
            ShowInTaskbar = false,
            ShowActivated = false,
            Visibility = Visibility.Hidden
        };
        window.Show();

        // Force layout
        parent.UpdateLayout();
        plotView.UpdateLayout();

        // Let dispatcher process layout/render
        System.Windows.Threading.Dispatcher.CurrentDispatcher.Invoke(
            System.Windows.Threading.DispatcherPriority.Background,
            new Action(delegate { }));

        var method = typeof(ChimeraAnalysisTabViewModel).GetMethod("CombinePlotAndLegend", BindingFlags.NonPublic | BindingFlags.Instance);
        Assert.That(method, Is.Not.Null);

        // Act
        var bitmap = method.Invoke(vm, new object[] { plotView, null }) as Bitmap;

        // Assert
        Assert.That(bitmap, Is.Not.Null);
        Assert.That(bitmap.Width, Is.EqualTo(1250));
        Assert.That(bitmap.Height, Is.EqualTo(625));
        bitmap.Dispose();
    }

    [Test]
    public void FindCommonAncestor_ReturnsNull_WhenNoCommonAncestor()
    {
        // Arrange
        var method = typeof(ChimeraAnalysisTabViewModel).GetMethod("FindCommonAncestor", BindingFlags.NonPublic | BindingFlags.Static);
        Assert.That(method, Is.Not.Null);

        var a = new Canvas();
        var b = new Canvas();

        // Act
        var ancestor = method.Invoke(null, new object[] { a, b });

        // Assert
        Assert.That(ancestor, Is.Null);
    }

    [Test]
    [NonParallelizable]
    public void ConstructChimericPsms_DoesNotCrash_WhenDeconDoesNotMatchPsms()
    {
        // Arrange: Copy a valid PSM and tweak its scan number values
        var allPsms = ChimeraGroupViewModelTests.AllMatchesMutable;
        var ms1ScanNumField = allPsms.First().GetType().GetProperty("PrecursorScanNum");
        foreach (var psm in allPsms)
        {
            // Use reflection to increment the public get, protected set PrecursorScanNum property by one
            var currentValue = (int)ms1ScanNumField.GetValue(psm);
            ms1ScanNumField.SetValue(psm, currentValue + 1, null);
        }

        var dataFiles = new Dictionary<string, MsDataFile>
        {
            { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile }
        };

        // Act: Should not throw even through many of the scans being wrong will throw exceptions
        List<ChimeraGroupViewModel> result = null;
        Assert.DoesNotThrow(() => result = ChimeraAnalysisTabViewModel.ConstructChimericPsms(allPsms, dataFiles));
    }

    [Test]
    [NonParallelizable]
    public void ConstructChimericPsms_SkipsGroup_WhenMs1OrMs2ScanIsMissing()
    {
        // Arrange: Copy a valid PSM and set its scan numbers to non-existent values
        var allPsms = ChimeraGroupViewModelTests.AllMatchesMutable;
        var ms1ScanNumField = allPsms.First().GetType().GetProperty("PrecursorScanNum");
        foreach (var psm in allPsms)
        {
            // Use reflection to increment the public get, protected set PrecursorScanNum property by one
            var currentValue = (int)ms1ScanNumField.GetValue(psm);
            ms1ScanNumField.SetValue(psm, currentValue + 184273, null);
        }

        var dataFiles = new Dictionary<string, MsDataFile>
        {
            { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile }
        };

        // Act: Should not throw even through many of the scans being wrong will throw exceptions
        List<ChimeraGroupViewModel> result = null;
        Assert.DoesNotThrow(() => result = ChimeraAnalysisTabViewModel.ConstructChimericPsms(allPsms, dataFiles));
        Assert.That(result, Is.Empty);
    }
}
