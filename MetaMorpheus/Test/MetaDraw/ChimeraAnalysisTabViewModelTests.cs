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

        var chimeraAnalysisTabViewModel = new ChimeraAnalysisTabViewModel(ChimeraGroupViewModelTests.AllMatches, dataFiles, TestExportDirectory);
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

        chimeraAnalysisTabViewModel = new ChimeraAnalysisTabViewModel(ChimeraGroupViewModelTests.AllMatches, dataFiles, TestExportDirectory);
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
    public static void ChimeraTabViewModel_PrecursorAssignmentIsCorrectInGroups()
    {
        var dataFiles = new Dictionary<string, MsDataFile>()
        {
            {"FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile }
        };
        var chimeraAnalysisTabViewModel = new ChimeraAnalysisTabViewModel(ChimeraGroupViewModelTests.AllMatches, dataFiles, TestExportDirectory);

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
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles, exportDir);

        // Assert
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
    public void ChimeraDrawnSequence_Setter_Notifies()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
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
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
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
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
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
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
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
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
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
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
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
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
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
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
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
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
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
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
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
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
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

    [Test]
    public void ConstructChimericPsms_ContinuesIfDataFileNotFound()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        // Provide an empty dataFiles dictionary so no data file can be found
        var dataFiles = new Dictionary<string, MsDataFile>();

        // Act
        var result = typeof(ChimeraAnalysisTabViewModel)
            .GetMethod("ConstructChimericPsms", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static)!
            .Invoke(null, new object[] { allPsms, dataFiles }) as List<ChimeraGroupViewModel>;

        // Assert
        Assert.That(result, Is.Empty, "No groups should be constructed if no data file is found.");
    }

    [Test]
    public void PrecursorIonAnnotations_AreLetterOnly_WhenUseLetterOnlyIsTrue()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles);

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
    public void ExportMs1_ThrowsOnUnsupportedExportType()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles);
        vm.SelectedChimeraGroup = vm.ChimeraGroupViewModels[0];
        vm.Ms1ChimeraPlot = new Ms1ChimeraPlot(new OxyPlot.Wpf.PlotView(), vm.SelectedChimeraGroup);
        vm.SelectedExportType = "UnsupportedType";

        // Act & Assert
        Assert.Throws<ArgumentOutOfRangeException>(() => vm.ExportMs1Command.Execute(null));
    }

    [Test]
    public void ExportMs2_ThrowsOnUnsupportedExportType()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles);
        vm.SelectedChimeraGroup = vm.ChimeraGroupViewModels[0];
        vm.ChimeraSpectrumMatchPlot = new ChimeraSpectrumMatchPlot(new OxyPlot.Wpf.PlotView(), vm.SelectedChimeraGroup);
        vm.SelectedExportType = "UnsupportedType";

        // Act & Assert
        Assert.Throws<ArgumentOutOfRangeException>(() => vm.ExportMs2Command.Execute(null));
    }
}
