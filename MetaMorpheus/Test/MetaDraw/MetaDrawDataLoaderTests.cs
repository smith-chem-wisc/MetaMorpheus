using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using GuiFunctions.MetaDraw;
using NUnit.Framework;
using Readers;
using GuiFunctions;
using MassSpectrometry;

namespace Test.MetaDraw;

[TestFixture]
public class MetaDrawDataLoaderTests
{
    private string mzmlPath;
    private MsDataFile msDataFile;
    private MetaDrawLogic logic;
    private MetaDrawDataLoader loader;

    [OneTimeSetUp]
    public void OneTimeSetup()
    {
        mzmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "SmallCalibratible_Yeast.mzML");

        Assert.That(File.Exists(mzmlPath), $"Test mzML file not found: {mzmlPath}");
        msDataFile = MsDataFileReader.GetDataFile(mzmlPath).LoadAllStaticData();

        logic = new MetaDrawLogic();
        loader = new MetaDrawDataLoader(logic);
    }

    [SetUp]
    public void SetUp()
    {
        logic.CleanUpResources();
    }

    [Test]
    public async Task LoadSpectraAsync_LoadsSpectraFiles()
    {
        logic.SpectraFilePaths.Add(mzmlPath);
        var errors = await loader.LoadSpectraAsync(CancellationToken.None);
        Assert.That(errors, Is.Empty);
        Assert.That(logic.MsDataFiles.Count, Is.EqualTo(1));
        Assert.That(logic.MsDataFiles.Values.First().GetMsDataScans().Count(), Is.GreaterThan(0));
    }

    [Test]
    public async Task LoadSpectraAsync_Cancellation_Throws()
    {
        var cts = new CancellationTokenSource();
        cts.Cancel();
        Assert.ThrowsAsync<TaskCanceledException>(async () => await loader.LoadSpectraAsync(cts.Token));
    }

    [Test]
    public async Task LoadPsmsAsync_NoFiles_ReturnsEmpty()
    {
        var errors = await loader.LoadPsmsAsync(false, CancellationToken.None);
        Assert.That(errors, Is.Empty);
    }

    [Test]
    public async Task LoadLibrariesAsync_NoFiles_ReturnsEmpty()
    {
        var errors = await loader.LoadLibrariesAsync(CancellationToken.None);
        Assert.That(errors, Is.Empty);
    }

    [Test]
    public async Task LoadAllAsync_SpectraOnly_LoadsSpectra()
    {
        logic.SpectraFilePaths.Add(mzmlPath);
        var errors = await loader.LoadAllAsync(true, false, false);
        Assert.That(errors, Is.Empty);
        Assert.That(logic.MsDataFiles.Count, Is.EqualTo(1));
    }

    [Test]
    public void ProcessBioPolymerTab_EnablesTabAndSetsExportDirectory()
    {
        var bioPolymerTab = new BioPolymerTabViewModel(logic, "C:\\Export");
        logic.AllSpectralMatches.Add(new DummySpectralmatch());
        loader.GetType().GetMethod("ProcessBioPolymerTab", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
            ?.Invoke(loader, new object[] { bioPolymerTab, "C:\\Export", CancellationToken.None });
        Assert.That(bioPolymerTab.IsTabEnabled, Is.True);
        Assert.That(bioPolymerTab.ExportDirectory, Is.EqualTo("C:\\Export"));
    }

    [Test]
    public void ProcessBioPolymerTab_LoadsDatabase()
    {
        var bioPolymerTab = new BioPolymerTabViewModel(logic, "C:\\Export");
        bioPolymerTab.DatabasePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
            @"TestData\TaGe_SA_A549_3_snip.fasta");
        logic.AllSpectralMatches.Add(new DummySpectralmatch());
        loader.GetType().GetMethod("ProcessBioPolymerTab", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
            ?.Invoke(loader, new object[] { bioPolymerTab, "C:\\Export", CancellationToken.None });
        Assert.That(bioPolymerTab.IsTabEnabled, Is.True);
        Assert.That(bioPolymerTab.ExportDirectory, Is.EqualTo("C:\\Export"));
    }

    [Test]
    public void ProcessDeconTab_EnablesTabAndSetsExportDirectory()
    {
        var deconTab = new DeconExplorationTabViewModel();
        logic.MsDataFiles.Clear();
        logic.MsDataFiles.Add("test", msDataFile);
        loader.GetType().GetMethod("ProcessDeconTab", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
            ?.Invoke(loader, new object[] { deconTab, "C:\\Export", CancellationToken.None });
        Assert.That(deconTab.IsTabEnabled, Is.True);
        Assert.That(deconTab.ExportDirectory, Is.EqualTo("C:\\Export"));
        Assert.That(deconTab.MsDataFiles.Count, Is.EqualTo(1));
    }

    [Test]
    public void ProcessChimeraTab_EnablesTabAndSetsExportDirectory()
    {
        var chimeraTab = new ChimeraAnalysisTabViewModel("C:\\Export");
        logic.FilteredListOfPsms.Add(new DummySpectralmatch());
        logic.MsDataFiles.Clear();
        logic.MsDataFiles.Add("test", msDataFile);
        loader.GetType().GetMethod("ProcessChimeraTab", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
            ?.Invoke(loader, new object[] { chimeraTab, "C:\\Export", CancellationToken.None });
        Assert.That(chimeraTab.IsTabEnabled, Is.True);
        Assert.That(chimeraTab.ExportDirectory, Is.EqualTo("C:\\Export"));
    }

    [Test]
    public void PrepareProgressVisualization_AddsStepsCorrectly()
    {
        logic.SpectraFilePaths.Add("file1.raw");
        logic.SpectraFilePaths.Add("file2.raw");
        logic.SpectralMatchResultFilePaths.Add("psm1.psmtsv");
        logic.SpectralLibraryPaths.Add("lib1.msp");
        loader.GetType().GetMethod("PrepareProgressVisualization", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
            ?.Invoke(loader, new object[] { true, true, true, logic });
        var vm = LoadingProgressViewModel.Instance;
        Assert.That(vm.Steps.Any(s => s.StepName == "Loading Spectra Files"), Is.True);
        Assert.That(vm.Steps.Any(s => s.StepName == "Loading Search Results"), Is.True);
        Assert.That(vm.Steps.Any(s => s.StepName == "Loading Spectral Libraries"), Is.True);
    }

    [Test]
    public void ReportProgress_UpdatesStep()
    {
        loader.GetType().GetMethod("AddStep", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
            ?.Invoke(loader, new object[] { "TestStep", 5 });
        loader.GetType().GetMethod("ReportProgress", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
            ?.Invoke(loader, new object[] { "TestStep", 3, 5 });
        var vm = LoadingProgressViewModel.Instance;
        var step = vm.Steps.FirstOrDefault(s => s.StepName == "TestStep");
        Assert.That(step, Is.Not.Null);
        Assert.That(step.Current, Is.EqualTo(3));
        Assert.That(step.Total, Is.EqualTo(5));
    }

    [Test]
    public void ProcessBioPolymerTab_Cancelled_DoesNotEnableTabOrProcess()
    {
        logic.SpectraFilePaths.Add("file1.raw");
        logic.SpectraFilePaths.Add("file2.raw");
        logic.SpectralMatchResultFilePaths.Add("psm1.psmtsv");
        logic.SpectralLibraryPaths.Add("lib1.msp");
        var tab = new BioPolymerTabViewModel(logic, "C:\\Export");
        tab.IsTabEnabled = false;
        tab.IsDatabaseLoaded = true;
        var cts = new CancellationTokenSource();
        cts.Cancel();
        loader.GetType().GetMethod("ProcessBioPolymerTab", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
            ?.Invoke(loader, new object[] { tab, "C:\\Export", cts.Token });
        Assert.That(tab.IsTabEnabled, Is.False);
    }

    [Test]
    public void ProcessDeconTab_Cancelled_DoesNotEnableTabOrProcess()
    {
        logic.SpectraFilePaths.Add("file1.raw");
        logic.SpectraFilePaths.Add("file2.raw");
        logic.SpectralMatchResultFilePaths.Add("psm1.psmtsv");
        logic.SpectralLibraryPaths.Add("lib1.msp");
        var tab = new DeconExplorationTabViewModel();
        tab.IsTabEnabled = false;
        var cts = new CancellationTokenSource();
        cts.Cancel();
        loader.GetType().GetMethod("ProcessDeconTab", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
            ?.Invoke(loader, new object[] { tab, "C:\\Export", cts.Token });
        Assert.That(tab.IsTabEnabled, Is.False);
    }

    [Test]
    public void ProcessChimeraTab_Cancelled_DoesNotEnableTabOrProcess()
    {
        logic.SpectraFilePaths.Add("file1.raw");
        logic.SpectraFilePaths.Add("file2.raw");
        logic.SpectralMatchResultFilePaths.Add("psm1.psmtsv");
        logic.SpectralLibraryPaths.Add("lib1.msp");
        var tab = new ChimeraAnalysisTabViewModel("C:\\Export");
        tab.IsTabEnabled = false;
        var cts = new CancellationTokenSource();
        cts.Cancel();
        loader.GetType().GetMethod("ProcessChimeraTab", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
            ?.Invoke(loader, new object[] { tab, "C:\\Export", cts.Token });
        Assert.That(tab.IsTabEnabled, Is.False);
    }

    [Test]
    [NonParallelizable]
    public void LoadAllAsync_FiresTabProcessingAndEnablesTabs()
    {
        if (!EverythingRunnerEngineTestCase.TryGetTestCase(EverythingRunnerEngineTestCases.BottomUpQValue, out var searchTestCase))
            Assert.Fail();

        string specFile1 = Path.Combine(TestContext.CurrentContext.TestDirectory,
            @"TestData\TaGe_SA_A549_3_snip.mzML");
        string specFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory,
            @"TestData\TaGe_SA_A549_3_snip_2.mzML");
        var dbPath = searchTestCase.DatabaseList.First().FilePath;
        var smPath = Directory.GetFiles(searchTestCase.OutputDirectory, "*.psmtsv", SearchOption.AllDirectories).First();

        logic.SpectraFilePaths.Add(specFile1);
        logic.SpectraFilePaths.Add(specFile2);
        logic.SpectralMatchResultFilePaths.Add(smPath);
        var bioTab = new BioPolymerTabViewModel(logic, "C:\\Export");
        bioTab.DatabasePath = dbPath;
        var deconTab = new DeconExplorationTabViewModel();
        var chimeraTab = new ChimeraAnalysisTabViewModel("C:\\Export");

        var errors = loader.LoadAllAsync(
            loadSpectra: true,
            loadPsms: true,
            loadLibraries: true,
            chimeraTabViewModel: chimeraTab,
            bioPolymerTabViewModel: bioTab,
            deconExplorationTabViewModel: deconTab);
        errors.Wait(); // wait for base processing to finish

        Thread.Sleep(5000); // wait for the asynchronous tab processing to finish

        Assert.That(errors.Result, Is.Empty);
        Assert.That(bioTab.IsTabEnabled, Is.True);
        Assert.That(deconTab.IsTabEnabled, Is.True);
        Assert.That(chimeraTab.IsTabEnabled, Is.True);
    }
}
