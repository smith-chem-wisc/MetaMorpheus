using System;
using System.IO;
using System.Linq;
using Avalonia;
using Avalonia.Headless;
using Avalonia.Headless.NUnit;
using EngineLayer;
using MetaMorpheus.Avalonia;
using MetaMorpheus.Avalonia.ViewModels;
using MetaMorpheus.Avalonia.Views;
using NUnit.Framework;

[assembly: AvaloniaTestApplication(typeof(Test.Avalonia.TestAppBuilder))]

namespace Test.Avalonia;

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
}

public class SpectraFileHandlingTests
{
    [TestCase("run.mzML", true)]
    [TestCase("run.mzml", true)]
    [TestCase("run.mgf", true)]
    [TestCase("run.raw", true)]
    [TestCase("notes.txt", false)]
    [TestCase("results.psmtsv", false)]
    public void OnlySpectraFormatsAreAccepted(string fileName, bool expected)
    {
        Assert.That(MainWindowViewModel.IsSupportedSpectraFile(fileName), Is.EqualTo(expected));
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
[Category("LongRunning")]
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

    /// <summary>RunCommand is async; wait for it rather than sleeping a fixed amount.</summary>
    private static void WaitForRunToFinish(MainWindowViewModel viewModel)
    {
        DateTime deadline = DateTime.UtcNow.AddMinutes(10);
        while (viewModel.IsRunning && DateTime.UtcNow < deadline)
        {
            System.Threading.Thread.Sleep(500);
        }

        Assert.That(viewModel.IsRunning, Is.False, "the search did not finish within 10 minutes");
    }
}
