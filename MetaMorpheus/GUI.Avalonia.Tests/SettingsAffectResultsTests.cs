using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using MetaMorpheusAvalonia.ViewModels;
using NUnit.Framework;
using TaskLayer;

namespace Test.AvaloniaGui;

/// <summary>
/// A settings dialog can look right, store values correctly, and still have no effect if the edited
/// task is not the one that runs. These drive real searches and compare the results, which is the only
/// way to show the dialog is wired all the way through to the engines.
/// </summary>
public class SettingsAffectResultsTests
{
    /// <summary>
    /// The default-settings search is the comparison point for every test here, so it runs once rather
    /// than per test. Three tests recomputing it meant six searches per platform leg where four do.
    /// </summary>
    private static int _baseline;

    private string _testData;
    private readonly List<string> _keptForDiagnosis = new();

    [SetUp]
    public void SetUp() => _testData = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");

    [OneTimeSetUp]
    public void RunTheBaselineOnce()
    {
        _testData = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
        _baseline = RunSearch(_ => { });
    }

    [TearDown]
    public void ReportKeptOutput()
    {
        foreach (string folder in _keptForDiagnosis)
        {
            TestContext.Out.WriteLine($"output kept for diagnosis: {folder}");
        }
        _keptForDiagnosis.Clear();
    }

    /// <summary>
    /// Runs a search with settings applied through the dialog's view model, returning the PSM count.
    /// Fails rather than returning a sentinel when no results are written: a sentinel of -1 satisfied
    /// every Is.LessThan assertion below, so a totally broken run read as a passing test.
    /// </summary>
    private int RunSearch(Action<TaskSettingsViewModel> configure)
    {
        var task = new SearchTask();
        var settings = new TaskSettingsViewModel(task, "Search");
        configure(settings);
        Assert.That(settings.Validate(), Is.Empty, "the test should only apply valid settings");
        settings.Apply();

        string outputFolder = Path.Combine(Path.GetTempPath(), "mm-settings-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(outputFolder);

        bool wrote = false;
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

            Assert.That(psmFile, Is.Not.Null,
                $"the search wrote no AllPSMs.psmtsv, so there is no result to compare. Output kept at {outputFolder}");

            wrote = true;
            return File.ReadAllLines(psmFile).Length - 1;
        }
        finally
        {
            // Keep the output when the run failed - deleting it unconditionally left nothing to
            // diagnose from a CI failure, which is exactly when it is needed.
            if (wrote && Directory.Exists(outputFolder))
            {
                Directory.Delete(outputFolder, recursive: true);
            }
            else if (!wrote)
            {
                _keptForDiagnosis.Add(outputFolder);
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
        Assert.That(_baseline, Is.EqualTo(22),
            "SearchTaskTest expects 22 target PSMs for mouseOne at default settings");
    }

    [Test]
    public void TighteningThePrecursorToleranceChangesTheResult()
    {
        int tightened = RunSearch(settings =>
        {
            settings.PrecursorTolerance = "0.05";
            settings.PrecursorToleranceIsPpm = true;
        });

        Assert.That(tightened, Is.Not.EqualTo(_baseline),
            $"a 0.05 ppm precursor tolerance should not give the same {_baseline} PSMs as the 5 ppm default; "
            + "if it does, the dialog's tolerance is not reaching the engine");
        Assert.That(tightened, Is.LessThan(_baseline), "a tighter tolerance should not find more PSMs");
    }

    [Test]
    public void RaisingTheMinimumPeptideLengthChangesTheResult()
    {
        int longPeptidesOnly = RunSearch(settings => settings.MinPeptideLength = 30);

        Assert.That(longPeptidesOnly, Is.LessThan(_baseline),
            "requiring 30-residue peptides should find fewer than the default minimum of 7");
    }

    /// <summary>
    /// A modification chosen in the dialog has to change what the search finds, not merely appear in
    /// CommonParameters. Fixed Carbamidomethyl is on by default, so removing it must move the count.
    /// </summary>
    [Test]
    public void RemovingADefaultFixedModificationChangesTheResult()
    {
        int withoutCarbamidomethyl = RunSearch(settings =>
        {
            foreach (ModificationChoice choice in settings.Modifications.Groups
                .SelectMany(g => g.Choices)
                .Where(c => c.IsFixed)
                .ToList())
            {
                choice.IsFixed = false;
            }
        });

        Assert.That(withoutCarbamidomethyl, Is.Not.EqualTo(_baseline),
            "dropping every fixed modification should change the result; if it does not, the "
            + "modification picker is not reaching the engine");
    }
}
