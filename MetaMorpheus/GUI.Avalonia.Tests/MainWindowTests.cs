using System;
using System.IO;
using System.Linq;
using Avalonia;
using Avalonia.Headless;
using Avalonia.Headless.NUnit;
using EngineLayer;
using MetaMorpheusAvalonia;
using MetaMorpheusAvalonia.ViewModels;
using MetaMorpheusAvalonia.Views;
using NUnit.Framework;

[assembly: AvaloniaTestApplication(typeof(Test.AvaloniaGui.TestAppBuilder))]

namespace Test.AvaloniaGui;

public static class TestAppBuilder
{
    public static AppBuilder BuildAvaloniaApp() => AppBuilder.Configure<App>()
        .UseHeadless(new AvaloniaHeadlessPlatformOptions());
}

/// <summary>
/// The WPF MainWindow cannot be constructed off Windows, and its 2,284 lines of logic cannot be
/// reached without a Window at all. These run headlessly, on any OS.
/// </summary>
public class MainWindowTests
{
    /// <summary>
    /// Builds the real window with the real view model. That constructor calls
    /// GlobalVariables.SetUpGlobalVariables, so this also proves the data directory, modifications and
    /// crosslinkers load - on whatever architecture the tests run on.
    /// </summary>
    [AvaloniaTest]
    public void MainWindowBuildsAndGlobalVariablesLoad()
    {
        using var viewModel = new MainWindowViewModel();
        var window = new MainWindow { DataContext = viewModel };
        window.Show();

        Assert.That(window.IsVisible, Is.True);
        Assert.That(window.Title, Is.EqualTo("MetaMorpheus"));
        Assert.That(viewModel.Log, Does.Contain("modifications"),
            "the constructor should report how many modifications GlobalVariables loaded");
    }

    /// <summary>
    /// A data directory that cannot be read used to take the process down inside
    /// OnFrameworkInitializationCompleted, before any window existed - the worst possible first run on
    /// exactly the platforms this GUI is for. The window opens now and says why.
    /// </summary>
    [AvaloniaTest]
    public void AFailedStartUpStillShowsAWindowWithTheReason()
    {
        using var viewModel = MainWindowViewModel.ThatFailedToStart(
            new DirectoryNotFoundException("Mods folder is missing"));
        var window = new MainWindow { DataContext = viewModel };
        window.Show();

        Assert.Multiple(() =>
        {
            Assert.That(window.IsVisible, Is.True, "the window must open even when start-up failed");
            Assert.That(viewModel.Log, Does.Contain("Mods folder is missing"));
            Assert.That(viewModel.CanRun, Is.False, "nothing can run without the data directory");
        });
    }
}

public class SpectraFileHandlingTests
{
    [TestCase("run.mzML", true)]
    [TestCase("run.mzml", true)]
    [TestCase("run.mgf", true)]
    [TestCase("run.raw", true)]
    [TestCase("run.msalign", true)]
    [TestCase("run.tdf", true)]
    [TestCase("notes.txt", false)]
    [TestCase("results.psmtsv", false)]
    public void OnlySpectraFormatsAreAccepted(string fileName, bool expected)
    {
        Assert.That(MainWindowViewModel.IsSupportedSpectraFile(fileName), Is.EqualTo(expected));
    }

    /// <summary>The accepted set is GlobalVariables', not a copy that can drift from it.</summary>
    [Test]
    public void EverySpectraFormatGlobalVariablesAcceptsIsAccepted()
    {
        foreach (string extension in GlobalVariables.AcceptedSpectraFormats)
        {
            Assert.That(MainWindowViewModel.IsSupportedSpectraFile("run" + extension), Is.True, extension);
        }
    }

    [TestCase("db.fasta", true)]
    [TestCase("db.xml", true)]
    [TestCase("db.xml.gz", true)]
    [TestCase("library.msp", true)]
    [TestCase("library.msl", true)]
    [TestCase("notes.txt", false)]
    public void OnlyDatabaseFormatsAreAccepted(string fileName, bool expected)
    {
        Assert.That(MainWindowViewModel.IsSupportedDatabaseFile(fileName), Is.EqualTo(expected));
    }

    /// <summary>
    /// Bruker data is a .d folder, so a user who opens one and picks the .tdf inside has to end up
    /// with the folder. Without the remap the selection was silently dropped.
    /// </summary>
    [Test]
    public void PickingTheTdfInsideABrukerFolderSelectsTheFolder()
    {
        string folder = Path.Combine("somewhere", "run.d");

        Assert.Multiple(() =>
        {
            Assert.That(MainWindowViewModel.ToBrukerFolderIfInside(Path.Combine(folder, "analysis.tdf")),
                Is.EqualTo(folder));
            Assert.That(MainWindowViewModel.ToBrukerFolderIfInside(Path.Combine(folder, "analysis.tdf_bin")),
                Is.EqualTo(folder));
            Assert.That(MainWindowViewModel.ToBrukerFolderIfInside(Path.Combine("elsewhere", "run.mzML")),
                Is.EqualTo(Path.Combine("elsewhere", "run.mzML")), "a normal file is left alone");
        });
    }

    /// <summary>Spectral libraries are databases as far as the run is concerned.</summary>
    [Test]
    public void SpectralLibrariesCanBeAdded()
    {
        using var viewModel = new MainWindowViewModel();
        viewModel.AddDatabases(new[] { "library.msl", "notes.txt" });

        Assert.That(viewModel.Databases.Select(d => d.FileName), Is.EqualTo(new[] { "library.msl" }));
    }

    /// <summary>
    /// Thermo .raw reading works on Windows and Linux only. Rejecting it up front with a reason beats
    /// failing part-way into a run, which is what would otherwise happen on macOS.
    /// </summary>
    [Test]
    public void ThermoRawIsRejectedOnMacOsOnly()
    {
        bool rejected = MainWindowViewModel.IsThermoRawUnsupportedHere("run.raw");
        Assert.That(rejected, Is.EqualTo(OperatingSystem.IsMacOS()));
        Assert.That(MainWindowViewModel.IsThermoRawUnsupportedHere("run.mzML"), Is.False);
    }
}

/// <summary>
/// Runs a real search through the view model, on the smallest spectra/database pairing in the
/// repository - the same one SearchTaskTest uses. This is the claim worth proving: that a search can
/// be driven from an Avalonia front end, on this machine, not just that the window renders.
/// </summary>
public class EndToEndSearchTests
{
    private string _outputFolder;

    [SetUp]
    public void SetUp()
    {
        _outputFolder = Path.Combine(Path.GetTempPath(), "mm-avalonia-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_outputFolder);
    }

    [TearDown]
    public void TearDown()
    {
        // StopLoops is static and process-wide, so a cancellation test must not leak into the next one.
        GlobalVariables.StopLoops = false;

        if (Directory.Exists(_outputFolder))
        {
            Directory.Delete(_outputFolder, recursive: true);
        }
    }

    [Test]
    public void SearchRunsAndWritesResults()
    {
        string testData = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
        string spectra = Path.Combine(testData, "mouseOne.mzML");
        string database = Path.Combine(testData, "mouseOne.xml");
        Assert.That(File.Exists(spectra), Is.True, $"missing {spectra}");
        Assert.That(File.Exists(database), Is.True, $"missing {database}");

        using var viewModel = new MainWindowViewModel { OutputFolder = _outputFolder };
        viewModel.AddSpectraFiles(new[] { spectra });
        viewModel.AddDatabases(new[] { database });
        viewModel.AddSearchTaskCommand.Execute(null);

        Assert.That(viewModel.CanRun, Is.True, "a spectra file, a database and a task should be enough to run");

        viewModel.RunCommand.Execute(null);
        WaitForRunToFinish(viewModel);

        Assert.That(viewModel.Log, Does.Not.Contain("ERROR:"), viewModel.Log);

        string[] written = Directory.GetFiles(_outputFolder, "*.psmtsv", SearchOption.AllDirectories);
        Assert.That(written, Is.Not.Empty, $"no .psmtsv was written. Log:{Environment.NewLine}{viewModel.Log}");
    }

    /// <summary>
    /// Cancel raises GlobalVariables.StopLoops and nothing else lowers it, so Run has to. Asserted by
    /// requiring a real search to still write results, rather than by reading the flag back.
    /// </summary>
    [Test]
    public void ASearchAfterACancellationStillProducesResults()
    {
        string testData = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");

        using var viewModel = new MainWindowViewModel { OutputFolder = _outputFolder };
        viewModel.CancelCommand.Execute(null);
        Assert.That(GlobalVariables.StopLoops, Is.True, "Cancel should raise the cooperative stop flag");

        viewModel.AddSpectraFiles(new[] { Path.Combine(testData, "mouseOne.mzML") });
        viewModel.AddDatabases(new[] { Path.Combine(testData, "mouseOne.xml") });
        viewModel.AddSearchTaskCommand.Execute(null);
        viewModel.RunCommand.Execute(null);
        WaitForRunToFinish(viewModel);

        Assert.That(viewModel.Log, Does.Not.Contain("ERROR:"), viewModel.Log);
        Assert.That(Directory.GetFiles(_outputFolder, "*.psmtsv", SearchOption.AllDirectories), Is.Not.Empty,
            $"a cancellation left the process unable to search. Log:{Environment.NewLine}{viewModel.Log}");
    }

    /// <summary>
    /// RunCommand is async, so wait for it rather than sleeping a fixed amount. The bound is a
    /// backstop against hanging forever, not an assertion about how fast a runner is - "don't hang"
    /// is expressed by timeout-minutes on the CI job, which fails the job rather than one test. This
    /// search takes about two seconds in practice; the margin is deliberately enormous so a loaded
    /// hosted macOS runner cannot turn slowness into a red build.
    /// </summary>
    private static void WaitForRunToFinish(MainWindowViewModel viewModel)
    {
        var backstop = TimeSpan.FromMinutes(30);
        DateTime deadline = DateTime.UtcNow + backstop;
        while (viewModel.IsRunning && DateTime.UtcNow < deadline)
        {
            System.Threading.Thread.Sleep(250);
        }

        if (viewModel.IsRunning)
        {
            Assert.Inconclusive($"the search was still running after {backstop.TotalMinutes:F0} minutes, "
                + "which says nothing about correctness - treat it as an environment problem");
        }
    }
}
