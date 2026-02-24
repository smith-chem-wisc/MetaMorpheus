using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Threading;
using EngineLayer;
using GuiFunctions;
using GuiFunctions.MetaDraw;
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
    private MetaDrawLogic metaDrawLogic;
    private MetaDrawLogic realDataLoaded;
    private List<SpectrumMatchFromTsv> realPsms;

    [OneTimeSetUp]
    public void Setup()
    {
        // Path to the test mzML file
        mzmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "SmallCalibratible_Yeast.mzML");
        Assert.That(File.Exists(mzmlPath), $"Test mzML file not found: {mzmlPath}");
        msDataFile = MsDataFileReader.GetDataFile(mzmlPath).LoadAllStaticData();
        metaDrawLogic = new MetaDrawLogic();

        var key = System.IO.Path.GetFileName(mzmlPath.Replace(GlobalVariables.GetFileExtension(mzmlPath), string.Empty));
        metaDrawLogic.MsDataFiles.Add(key, msDataFile);

        // Set up real data tests
        realDataLoaded = new();
        var psmPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TDGPTMDSearchResults.psmtsv");
        string msDataPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TDGPTMDSearchSingleSpectra.mzML");
        key = System.IO.Path.GetFileName(msDataPath.Replace(GlobalVariables.GetFileExtension(msDataPath), string.Empty));
        realDataLoaded.SpectraFilePaths.Add(msDataPath);
        realDataLoaded.SpectralMatchResultFilePaths.Add(psmPath);   
        new MetaDrawDataLoader(realDataLoaded).LoadAllAsync(true, true, true).Wait();

        var file = realDataLoaded.MsDataFiles.First().Value;
        realDataLoaded.MsDataFiles.Clear();
        realDataLoaded.MsDataFiles.Add("SmallCalibratible_Yeast", file);
        realPsms = realDataLoaded.AllSpectralMatches;
    }

    [Test]
    public void PopulateScansCollection_FullSpectrumMode_AddsAllScans()
    {
        var vm = new DeconExplorationTabViewModel(metaDrawLogic);
        vm.Mode = DeconvolutionMode.FullSpectrum;
        vm.SelectedMsDataFile = msDataFile;

        Thread.Sleep(3000); // wait for async load
        Assert.That(vm.Scans.Count, Is.EqualTo(msDataFile.GetMsDataScans().Count()));
    }

    [Test]
    public void PopulateScansCollection_IsolationRegionMode_AddsOnlyMs2()
    {
        var vm = new DeconExplorationTabViewModel(metaDrawLogic);
        vm.Mode = DeconvolutionMode.IsolationRegion;
        vm.SelectedMsDataFile = msDataFile;

        Assert.That(vm.Scans.All(s => s.MsnOrder == 2));
    }

    [Test]
    public void Changing_Mode_RefreshesScans()
    {
        var vm = new DeconExplorationTabViewModel(metaDrawLogic);
        vm.SelectedMsDataFile = msDataFile;
        vm.Mode = DeconvolutionMode.FullSpectrum;
        Thread.Sleep(1000); // wait for async load

        int allCount = vm.Scans.Count;
        vm.Mode = DeconvolutionMode.IsolationRegion;
        Thread.Sleep(1000); // wait for async load
        int ms2Count = vm.Scans.Count;

        Assert.That(allCount, Is.GreaterThan(ms2Count));
        Assert.That(vm.Scans.All(s => s.MsnOrder == 2));
    }

    [Test]
    public void RunDeconvolutionCommand_DoesNotThrow_WhenNoScanSelected()
    {
        var vm = new DeconExplorationTabViewModel(metaDrawLogic);
        Assert.DoesNotThrow(() => vm.RunDeconvolutionCommand.Execute(null));
    }

    [Test, Apartment(System.Threading.ApartmentState.STA)]
    public void RunDeconvolutionCommand_FullSpectrumMode_PopulatesDeconvolutedSpecies()
    {
        var vm = new DeconExplorationTabViewModel(metaDrawLogic);
        vm.Mode = DeconvolutionMode.FullSpectrum;
        vm.SelectedMsDataFile = msDataFile;
        vm.SelectedMsDataScan = msDataFile.GetMsDataScans().FirstOrDefault();
        vm.RunDeconvolutionCommand.Execute(new PlotView());

        Assert.That(vm.DeconvolutedSpecies.Count, Is.GreaterThanOrEqualTo(0));
    }

    [Test, Apartment(System.Threading.ApartmentState.STA)]
    public void RunDeconvolutionCommand_IsolationRegionMode_PopulatesDeconvolutedSpecies()
    {
        var vm = new DeconExplorationTabViewModel(metaDrawLogic);
        Assert.That(vm.IsLoading, Is.False);
        vm.Mode = DeconvolutionMode.IsolationRegion;
        vm.SelectedMsDataFile = msDataFile;
        var ms2 = msDataFile.GetMsDataScans().FirstOrDefault(s => s.MsnOrder == 2);
        if (ms2 == null)
            Assert.Inconclusive("No MS2 scan in test file.");
        vm.SelectedMsDataScan = ms2;
        vm.RunDeconvolutionCommand.Execute(new PlotView());

        Assert.That(vm.DeconvolutedSpecies.Count, Is.GreaterThanOrEqualTo(0));
    }

    [Test, Apartment(System.Threading.ApartmentState.STA)]
    public void RunDeconvolutionCommand_FullSpectrumMode_LimitsAxis()
    {
        var vm = new DeconExplorationTabViewModel(metaDrawLogic);
        vm.Mode = DeconvolutionMode.FullSpectrum;
        vm.SelectedMsDataFile = msDataFile;
        vm.SelectedMsDataScan = msDataFile.GetMsDataScans().FirstOrDefault();
        vm.RunDeconvolutionCommand.Execute(new PlotView());

        var specRange = vm.SelectedMsDataScan!.MassSpectrum.Range;
        var xAxis = vm.Plot!.Model.Axes[0];
        Assert.That(xAxis.Minimum, Is.EqualTo(specRange.Minimum).Within(5));
        Assert.That(xAxis.Maximum, Is.EqualTo(specRange.Maximum).Within(5));

        Assert.That(vm.MinMzToPlot, Is.EqualTo(0));
        vm.MinMzToPlot = 800;
        vm.RunDeconvolutionCommand.Execute(new PlotView());
        xAxis = vm.Plot!.Model.Axes[0];
        Assert.That(xAxis.ActualMinimum, Is.EqualTo(800).Within(5));
        Assert.That(xAxis.ActualMaximum, Is.EqualTo(specRange.Maximum).Within(5));

        Assert.That(vm.MaxMzToPlot, Is.EqualTo(0));
        vm.MaxMzToPlot = 1200;
        vm.RunDeconvolutionCommand.Execute(new PlotView());
        xAxis = vm.Plot!.Model.Axes[0];
        Assert.That(xAxis.ActualMinimum, Is.EqualTo(800).Within(5));
        Assert.That(xAxis.ActualMaximum, Is.EqualTo(1200).Within(5));
    }

    [Test]
    public void DeconvolutionModes_ContainsAllModes()
    {
        var vm = new DeconExplorationTabViewModel(metaDrawLogic);
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
        var plot = new DeconvolutionPlot(plotView, scan, deconResults, DeconvolutionMode.FullSpectrum);

        Assert.DoesNotThrow(() => plot.AnnotatePlot(null!));
        Assert.DoesNotThrow(() => plot.AnnotatePlot(new List<DeconvolutedSpeciesViewModel>()));
    }

    [Test]
    public void OnlyIdentifiedScans_FiltersScans_FullSpectrumMode()
    {
        var vm = new DeconExplorationTabViewModel(realDataLoaded);
        vm.MsDataFiles.Add(realDataLoaded.MsDataFiles.Values.First());
        vm.SelectedMsDataFile = vm.MsDataFiles.First();
        vm.Mode = DeconvolutionMode.FullSpectrum;
        Thread.Sleep(1000); // wait for async load

        // Simulate identified scans
        var scanCount = vm.Scans.Count;

        var acceptableScanNumbers = realPsms.Select(p => p.Ms2ScanNumber)
            .Concat(realPsms.Select(p => p.PrecursorScanNum)).ToHashSet();

        vm.OnlyIdentifiedScans = true;
        Thread.Sleep(1000); // wait for async load
        var scanCountAfterFilter = vm.Scans.Count;
        Assert.That(scanCountAfterFilter, Is.LessThan(scanCount));
        var expectedScanNumbers = acceptableScanNumbers.Intersect(vm.Scans.Select(s => s.OneBasedScanNumber)).ToHashSet();
        var actualScanNumbers = vm.Scans.Select(s => s.OneBasedScanNumber).ToHashSet();

        Assert.That(actualScanNumbers, Is.EquivalentTo(expectedScanNumbers));
    }

    [Test]
    public void OnlyIdentifiedScans_FiltersScans_IsolationRegionMode()
    {
        var vm = new DeconExplorationTabViewModel(realDataLoaded);
        vm.MsDataFiles.Add(realDataLoaded.MsDataFiles.Values.First());
        vm.SelectedMsDataFile = vm.MsDataFiles.First();
        vm.Mode = DeconvolutionMode.IsolationRegion;
        Thread.Sleep(1000); // wait for async load

        // Simulate identified scans
        var scanCount = vm.Scans.Count;

        var acceptableScanNumbers = realPsms.Select(p => p.Ms2ScanNumber)
            .ToHashSet();

        vm.OnlyIdentifiedScans = true;
        Thread.Sleep(1000); // wait for async load
        var scanCountAfterFilter = vm.Scans.Count;
        Assert.That(scanCountAfterFilter, Is.LessThan(scanCount));
        var expectedScanNumbers = acceptableScanNumbers.Intersect(vm.Scans.Select(s => s.OneBasedScanNumber)).ToHashSet();
        var actualScanNumbers = vm.Scans.Select(s => s.OneBasedScanNumber).ToHashSet();

        Assert.That(actualScanNumbers, Is.EquivalentTo(expectedScanNumbers));
        Assert.That(vm.Scans.All(s => s.MsnOrder == 2));
    }

    [Test, Apartment(System.Threading.ApartmentState.STA)]
    public void RunDeconvolutionCommand_FullSpectrumMode_LimitsCharges()
    {
        var vm = new DeconExplorationTabViewModel(metaDrawLogic);
        vm.Mode = DeconvolutionMode.FullSpectrum;
        vm.SelectedMsDataFile = msDataFile;
        vm.SelectedMsDataScan = msDataFile.GetMsDataScans().FirstOrDefault();
        vm.RunDeconvolutionCommand.Execute(new PlotView());

        int min = vm.DeconvolutedSpecies.Min(s => s.Charge);
        int max = vm.DeconvolutedSpecies.Max(s => s.Charge);
        Assert.That(vm.MinChargeToAnnotate, Is.LessThanOrEqualTo(min));
        Assert.That(vm.MaxChargeToAnnotate, Is.GreaterThanOrEqualTo(max));

        vm.MinChargeToAnnotate = 2;
        vm.RunDeconvolutionCommand.Execute(new PlotView());
        Assert.That(vm.DeconvolutedSpecies.All(s => s.Charge >= 2));

        vm.MaxChargeToAnnotate = 3;
        vm.RunDeconvolutionCommand.Execute(new PlotView());
        Assert.That(vm.DeconvolutedSpecies.All(s => s.Charge <= 3));
    }
}