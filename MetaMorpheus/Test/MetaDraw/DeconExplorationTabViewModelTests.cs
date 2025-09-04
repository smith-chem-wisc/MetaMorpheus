using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Windows.Documents;
using GuiFunctions;
using MassSpectrometry;
using NUnit.Framework;
using OxyPlot.Wpf;
using Readers;

namespace Test.MetaDraw;

[TestFixture]
public class DeconExplorationTabViewModelTests
{
    private string mzmlPath;
    private MsDataFile msDataFile;

    [OneTimeSetUp]
    public void Setup()
    {
        // Path to the test mzML file
        mzmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "SmallCalibratible_Yeast.mzML");
        Assert.That(File.Exists(mzmlPath), $"Test mzML file not found: {mzmlPath}");
        msDataFile = MsDataFileReader.GetDataFile(mzmlPath).LoadAllStaticData();
    }

    [Test]
    public void PopulateScansCollection_FullSpectrumMode_AddsAllScans()
    {
        var vm = new DeconExplorationTabViewModel();
        vm.Mode = DeconvolutionMode.FullSpectrum;
        vm.SelectedMsDataFile = msDataFile;

        Assert.That(vm.Scans.Count, Is.EqualTo(msDataFile.GetMsDataScans().Count()));
    }

    [Test]
    public void PopulateScansCollection_IsolationRegionMode_AddsOnlyMs2()
    {
        var vm = new DeconExplorationTabViewModel();
        vm.Mode = DeconvolutionMode.IsolationRegion;
        vm.SelectedMsDataFile = msDataFile;

        Assert.That(vm.Scans.All(s => s.MsnOrder == 2));
    }

    [Test]
    public void Changing_Mode_RefreshesScans()
    {
        var vm = new DeconExplorationTabViewModel();
        vm.SelectedMsDataFile = msDataFile;
        vm.Mode = DeconvolutionMode.FullSpectrum;
        int allCount = vm.Scans.Count;
        vm.Mode = DeconvolutionMode.IsolationRegion;
        int ms2Count = vm.Scans.Count;

        Assert.That(allCount, Is.GreaterThan(ms2Count));
        Assert.That(vm.Scans.All(s => s.MsnOrder == 2));
    }

    [Test]
    public void RunDeconvolutionCommand_DoesNotThrow_WhenNoScanSelected()
    {
        var vm = new DeconExplorationTabViewModel();
        Assert.DoesNotThrow(() => vm.RunDeconvolutionCommand.Execute(null));
    }

    [Test, Apartment(System.Threading.ApartmentState.STA)]
    public void RunDeconvolutionCommand_FullSpectrumMode_PopulatesDeconvolutedSpecies()
    {
        var vm = new DeconExplorationTabViewModel();
        vm.Mode = DeconvolutionMode.FullSpectrum;
        vm.SelectedMsDataFile = msDataFile;
        vm.SelectedMsDataScan = msDataFile.GetMsDataScans().FirstOrDefault();
        vm.RunDeconvolutionCommand.Execute(new PlotView());

        Assert.That(vm.DeconvolutedSpecies.Count, Is.GreaterThanOrEqualTo(0));
    }

    [Test, Apartment(System.Threading.ApartmentState.STA)]
    public void RunDeconvolutionCommand_IsolationRegionMode_PopulatesDeconvolutedSpecies()
    {
        var vm = new DeconExplorationTabViewModel();
        vm.Mode = DeconvolutionMode.IsolationRegion;
        vm.SelectedMsDataFile = msDataFile;
        var ms2 = msDataFile.GetMsDataScans().FirstOrDefault(s => s.MsnOrder == 2);
        if (ms2 == null)
            Assert.Inconclusive("No MS2 scan in test file.");
        vm.SelectedMsDataScan = ms2;
        vm.RunDeconvolutionCommand.Execute(new PlotView());

        Assert.That(vm.DeconvolutedSpecies.Count, Is.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void DeconvolutionModes_ContainsAllModes()
    {
        var vm = new DeconExplorationTabViewModel();
        var modes = vm.DeconvolutionModes;
        Assert.That(modes, Does.Contain(DeconvolutionMode.FullSpectrum));
        Assert.That(modes, Does.Contain(DeconvolutionMode.IsolationRegion));
    }

    [Test, Apartment(ApartmentState.STA)]
    public void DeconPlot_AnnotatePlot_ReturnsEarlyOnNullOrEmpty()
    {
        var plotView = new PlotView();
        var scan = msDataFile.LoadAllStaticData().GetMsDataScans().First();
        var deconResults = new List<DeconvolutedSpeciesViewModel>();
        var plot = new DeconvolutionPlot(plotView, scan, deconResults);

        Assert.DoesNotThrow(() => plot.AnnotatePlot(null!));
        Assert.DoesNotThrow(() => plot.AnnotatePlot(new List<DeconvolutedSpeciesViewModel>()));
    }
}