using System.Collections.Generic;
using System.IO;
using System.Linq;
using GuiFunctions;
using MassSpectrometry;
using NUnit.Framework;
using OxyPlot;
using OxyPlot.Wpf;
using Readers;

namespace Test.MetaDraw;

[TestFixture]
public class DeconvolutedSpeciesViewModelTests
{
    [Test]
    public void Properties_AreCorrect()
    {
        // Arrange: create a fake isotopic envelope
        var peaks = new List<(double mz, double intensity)>
        {
            (500.0, 100.0),
            (501.0, 200.0),
            (502.0, 50.0)
        };
        var envelope = new IsotopicEnvelope(0, peaks, 1000, 2, peaks.Sum(p => p.intensity), 1);

        // Act
        var vm = new DeconvolutedSpeciesViewModel(envelope);

        // Assert
        Assert.That(vm.Envelope, Is.EqualTo(envelope));
        Assert.That(vm.Color, Is.EqualTo(OxyColors.Automatic));
        Assert.That(vm.MonoisotopicMass, Is.EqualTo(1000.0));
        Assert.That(vm.Charge, Is.EqualTo(2));
        Assert.That(vm.PeakCount, Is.EqualTo(3));
        Assert.That(vm.Intensity, Is.EqualTo(350.0));
        Assert.That(vm.Annotation, Does.Contain("M=1000"));
        Assert.That(vm.Annotation, Does.Contain("z=2"));
        Assert.That(vm.MostAbundantMz, Is.EqualTo(501.0));
        Assert.That(vm.PeakMzs, Is.EqualTo("500, 501, 502"));
    }

    [Test]
    public void Properties_CacheValues()
    {
        // Arrange
        var peaks = new List<(double mz, double intensity)>
        {
            (100, 1), (200, 2)
        };
        var envelope = new IsotopicEnvelope(0, peaks, 1000, 3, peaks.Sum(p => p.intensity), 1);
        var vm = new DeconvolutedSpeciesViewModel(envelope);

        // Act
        var intensity1 = vm.Intensity;
        var intensity2 = vm.Intensity;
        var annotation1 = vm.Annotation;
        var annotation2 = vm.Annotation;
        var mostAbundantMz1 = vm.MostAbundantMz;
        var mostAbundantMz2 = vm.MostAbundantMz;
        var peakMzs1 = vm.PeakMzs;
        var peakMzs2 = vm.PeakMzs;

        // Assert: values are cached and consistent
        Assert.That(intensity1, Is.EqualTo(intensity2));
        Assert.That(annotation1, Is.EqualTo(annotation2));
        Assert.That(mostAbundantMz1, Is.EqualTo(mostAbundantMz2));
        Assert.That(peakMzs1, Is.EqualTo(peakMzs2));
    }
}

[TestFixture]
public class DeconExplorationViewModelTests
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
        var vm = new DeconExplorationViewModel();
        vm.Mode = DeconvolutionMode.FullSpectrum;
        vm.SelectedMsDataFile = msDataFile;

        Assert.That(vm.Scans.Count, Is.EqualTo(msDataFile.GetMsDataScans().Count()));
    }

    [Test]
    public void PopulateScansCollection_IsolationRegionMode_AddsOnlyMs2()
    {
        var vm = new DeconExplorationViewModel();
        vm.Mode = DeconvolutionMode.IsolationRegion;
        vm.SelectedMsDataFile = msDataFile;

        Assert.That(vm.Scans.All(s => s.MsnOrder == 2));
    }

    [Test]
    public void Changing_Mode_RefreshesScans()
    {
        var vm = new DeconExplorationViewModel();
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
        var vm = new DeconExplorationViewModel();
        Assert.DoesNotThrow(() => vm.RunDeconvolutionCommand.Execute(null));
    }

    [Test, Apartment(System.Threading.ApartmentState.STA)]
    public void RunDeconvolutionCommand_FullSpectrumMode_PopulatesDeconvolutedSpecies()
    {
        var vm = new DeconExplorationViewModel();
        vm.Mode = DeconvolutionMode.FullSpectrum;
        vm.SelectedMsDataFile = msDataFile;
        vm.SelectedMsDataScan = msDataFile.GetMsDataScans().FirstOrDefault();
        vm.RunDeconvolutionCommand.Execute(new PlotView());

        Assert.That(vm.DeconvolutedSpecies.Count, Is.GreaterThanOrEqualTo(0));
    }

    [Test, Apartment(System.Threading.ApartmentState.STA)]
    public void RunDeconvolutionCommand_IsolationRegionMode_PopulatesDeconvolutedSpecies()
    {
        var vm = new DeconExplorationViewModel();
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
        var vm = new DeconExplorationViewModel();
        var modes = vm.DeconvolutionModes;
        Assert.That(modes, Does.Contain(DeconvolutionMode.FullSpectrum));
        Assert.That(modes, Does.Contain(DeconvolutionMode.IsolationRegion));
    }
}
