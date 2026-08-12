using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using MetaMorpheus.Avalonia.ViewModels;
using NUnit.Framework;
using TaskLayer;

namespace Test.Avalonia;

/// <summary>
/// A settings dialog can look right, store values correctly, and still have no effect if the edited
/// task is not the one that runs. These drive real searches and compare the results, which is the only
/// way to show the dialog is wired all the way through to the engines.
/// </summary>
[Category("LongRunning")]
public class SettingsAffectResultsTests
{
    private string _testData;

    [SetUp]
    public void SetUp()
    {
        GlobalVariables.SetUpGlobalVariables();
        _testData = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
    }

    /// <summary>Runs a search with settings applied through the dialog's view model, returning the PSM count.</summary>
    private int RunSearch(Action<TaskSettingsViewModel> configure)
    {
        var task = new SearchTask();
        var settings = new TaskSettingsViewModel(task, "Search");
        configure(settings);
        Assert.That(settings.Validate(), Is.Empty, "the test should only apply valid settings");
        settings.Apply();

        string outputFolder = Path.Combine(Path.GetTempPath(), "mm-settings-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(outputFolder);
        try
        {
            new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("Task1-Search", task) },
                new List<string> { Path.Combine(_testData, "mouseOne.mzML") },
                new List<DbForTask> { new DbForTask(Path.Combine(_testData, "mouseOne.xml"), false) },
                outputFolder).Run();

            string psmFile = Directory
                .GetFiles(outputFolder, "AllPSMs.psmtsv", SearchOption.AllDirectories)
                .FirstOrDefault();
            return psmFile is null ? -1 : File.ReadAllLines(psmFile).Length - 1;
        }
        finally
        {
            if (Directory.Exists(outputFolder))
            {
                Directory.Delete(outputFolder, recursive: true);
            }
        }
    }

    /// <summary>
    /// The default settings must reproduce the count MetaMorpheus's own SearchTaskTest asserts for this
    /// file, so a change in this GUI cannot quietly diverge from the established result.
    /// </summary>
    [Test]
    public void DefaultSettingsReproduceTheKnownPsmCount()
    {
        Assert.That(RunSearch(_ => { }), Is.EqualTo(22),
            "SearchTaskTest expects 22 target PSMs for mouseOne at default settings");
    }

    [Test]
    public void TighteningThePrecursorToleranceChangesTheResult()
    {
        int baseline = RunSearch(_ => { });
        int tightened = RunSearch(settings =>
        {
            settings.PrecursorTolerance = "0.05";
            settings.PrecursorToleranceIsPpm = true;
        });

        Assert.That(tightened, Is.Not.EqualTo(baseline),
            $"a 0.05 ppm precursor tolerance should not give the same {baseline} PSMs as the 5 ppm default; "
            + "if it does, the dialog's tolerance is not reaching the engine");
        Assert.That(tightened, Is.LessThan(baseline), "a tighter tolerance should not find more PSMs");
    }

    [Test]
    public void RaisingTheMinimumPeptideLengthChangesTheResult()
    {
        int baseline = RunSearch(_ => { });
        int longPeptidesOnly = RunSearch(settings => settings.MinPeptideLength = 30);

        Assert.That(longPeptidesOnly, Is.LessThan(baseline),
            "requiring 30-residue peptides should find fewer than the default minimum of 7");
    }
}
