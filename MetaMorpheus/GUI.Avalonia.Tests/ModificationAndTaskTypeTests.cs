using System.Linq;
using EngineLayer;
using EngineLayer.GlycoSearch;
using MetaMorpheusAvalonia.ViewModels;
using NUnit.Framework;
using TaskLayer;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;

namespace Test.AvaloniaGui;

/// <summary>
/// Modification selection is the setting most likely to be silently lost, because Apply() rebuilds
/// CommonParameters and the mods are carried as (ModificationType, IdWithMotif) pairs rather than
/// objects. These check the pairs survive the round trip.
/// </summary>
public class ModificationSelectionTests
{
    [Test]
    public void TheDefaultsAreLoadedAsSelected()
    {
        var settings = new TaskSettingsViewModel(new SearchTask(), "Search");

        Assert.That(settings.Modifications.TotalCount, Is.GreaterThan(1000),
            "GlobalVariables loads thousands of modifications, so all of them should be offered");
        Assert.That(settings.Modifications.FixedSelection.Select(m => m.Item2),
            Does.Contain("Carbamidomethyl on C"), "the default fixed modification should be ticked");
        Assert.That(settings.Modifications.VariableSelection.Select(m => m.Item2),
            Does.Contain("Oxidation on M"), "the default variable modification should be ticked");
    }

    [Test]
    public void SelectingAModificationReachesTheTask()
    {
        var task = new SearchTask();
        var settings = new TaskSettingsViewModel(task, "Search");

        ModificationChoice acetyl = settings.Modifications.Groups
            .SelectMany(g => g.Choices)
            .First(c => c.Name.StartsWith("Acetylation on K"));
        settings.Modifications.SetVariable(acetyl, true);

        settings.Apply();

        Assert.That(task.CommonParameters.ListOfModsVariable.Select(m => m.Item2),
            Does.Contain(acetyl.Name), "a newly ticked variable modification must reach the task");
    }

    [Test]
    public void DeselectingAModificationRemovesItFromTheTask()
    {
        var task = new SearchTask();
        var settings = new TaskSettingsViewModel(task, "Search");

        ModificationChoice oxidation = settings.Modifications.Groups
            .SelectMany(g => g.Choices)
            .First(c => c.Name == "Oxidation on M");
        settings.Modifications.SetVariable(oxidation, false);

        settings.Apply();

        Assert.That(task.CommonParameters.ListOfModsVariable.Select(m => m.Item2),
            Does.Not.Contain("Oxidation on M"), "unticking must remove it, not just fail to add it");
    }

    /// <summary>
    /// A modification cannot be both fixed and variable without the search treating that residue
    /// inconsistently, so choosing one clears the other.
    /// </summary>
    [Test]
    public void FixedAndVariableAreMutuallyExclusive()
    {
        var settings = new TaskSettingsViewModel(new SearchTask(), "Search");
        ModificationChoice choice = settings.Modifications.Groups.SelectMany(g => g.Choices).First();

        settings.Modifications.SetVariable(choice, true);
        settings.Modifications.SetFixed(choice, true);

        Assert.That(choice.IsFixed, Is.True);
        Assert.That(choice.IsVariable, Is.False, "marking it fixed should clear the variable tick");
    }

    /// <summary>
    /// With thousands of modifications the list has to be filterable, but a filter must never hide
    /// something already ticked, or the user loses a selection without seeing it happen.
    /// </summary>
    [Test]
    public void FilteringNeverHidesASelection()
    {
        var settings = new TaskSettingsViewModel(new SearchTask(), "Search");

        settings.Modifications.Filter = "zzzz-matches-nothing";

        var visible = settings.Modifications.Groups.SelectMany(g => g.Choices).ToList();
        Assert.That(visible, Is.Not.Empty, "the ticked defaults should still be shown");
        Assert.That(visible.All(c => c.IsFixed || c.IsVariable), Is.True,
            "with a filter matching nothing, only the selected modifications should remain visible");
    }

    [Test]
    public void FilteringNarrowsTheList()
    {
        var settings = new TaskSettingsViewModel(new SearchTask(), "Search");
        int unfiltered = settings.Modifications.Groups.SelectMany(g => g.Choices).Count();

        settings.Modifications.Filter = "Phospho";
        int filtered = settings.Modifications.Groups.SelectMany(g => g.Choices).Count();

        Assert.That(filtered, Is.LessThan(unfiltered));
        Assert.That(settings.Modifications.Groups.SelectMany(g => g.Choices)
            .Any(c => c.Name.Contains("Phospho")), Is.True);
    }

    /// <summary>
    /// The exclusivity rule has to hold however the property is set, not only through SetFixed. A
    /// CheckBox bound straight to IsFixed - the obvious binding - used to bypass it and produce the
    /// both-fixed-and-variable state the rule exists to prevent.
    /// </summary>
    [Test]
    public void ExclusivityHoldsWhenThePropertyIsSetDirectly()
    {
        var settings = new TaskSettingsViewModel(new SearchTask(), "Search");
        ModificationChoice choice = settings.Modifications.Groups.SelectMany(g => g.Choices).First();

        choice.IsVariable = true;
        choice.IsFixed = true;

        Assert.Multiple(() =>
        {
            Assert.That(choice.IsFixed, Is.True);
            Assert.That(choice.IsVariable, Is.False, "setting IsFixed directly must still clear IsVariable");
        });

        choice.IsVariable = true;
        Assert.That(choice.IsFixed, Is.False, "and the reverse");
    }

    /// <summary>Setting the property directly must also refresh the counts the header shows.</summary>
    [Test]
    public void CountsFollowADirectPropertySet()
    {
        var settings = new TaskSettingsViewModel(new SearchTask(), "Search");
        int before = settings.Modifications.SelectedVariableCount;

        settings.Modifications.Groups.SelectMany(g => g.Choices)
            .First(c => !c.IsFixed && !c.IsVariable).IsVariable = true;

        Assert.That(settings.Modifications.SelectedVariableCount, Is.EqualTo(before + 1));
    }

    /// <summary>
    /// Collapsing to one entry per (ModificationType, IdWithMotif) is forced - that pair is the
    /// identity TOML stores, so two modifications sharing it cannot be selected apart. It is only safe
    /// while the entries sharing a pair are interchangeable. On the shipped databases 84 pairs collapse,
    /// all in the glycan sets, and every one is a duplicate rather than a distinct modification. This
    /// fails if a database change ever makes that untrue, rather than quietly hiding a modification.
    /// </summary>
    [Test]
    public void CollapsingDuplicatesNeverHidesADistinctModification()
    {
        var selection = new TaskSettingsViewModel(new SearchTask(), "Search").Modifications;

        Assert.That(selection.DistinctModificationsCollapsed, Is.Zero,
            $"{selection.DistinctModificationsCollapsed} of {selection.CollapsedDuplicates} collapsed "
            + "identities hold modifications that differ in location, mass or target, so collapsing "
            + "them makes one of them unreachable");
        Assert.That(selection.CollapsedDuplicates, Is.GreaterThan(0),
            "the shipped databases do contain duplicates, so this is measuring something");
    }

    /// <summary>Defaults come from CommonParameters, not a copied literal, so they cannot drift.</summary>
    [Test]
    public void ResetToDefaultsMatchesTheEnginesOwnDefaults()
    {
        var settings = new TaskSettingsViewModel(new SearchTask(), "Search");
        ModificationChoice spare = settings.Modifications.Groups.SelectMany(g => g.Choices)
            .First(c => !c.IsFixed && !c.IsVariable);
        spare.IsVariable = true;

        settings.Modifications.ResetToDefaults();

        var engineDefaults = new CommonParameters();
        Assert.Multiple(() =>
        {
            Assert.That(settings.Modifications.FixedSelection, Is.EquivalentTo(engineDefaults.ListOfModsFixed));
            Assert.That(settings.Modifications.VariableSelection, Is.EquivalentTo(engineDefaults.ListOfModsVariable));
        });
    }

    /// <summary>
    /// CommonParameters picks a different default set for RNA, which the old hardcoded list could not
    /// express - it only ever produced the peptide defaults.
    /// </summary>
    [Test]
    public void ResetToDefaultsUsesTheRnaDefaultsForAnRnaTask()
    {
        var rnaDigestion = new RnaDigestionParams();
        var settings = new TaskSettingsViewModel(
            new SearchTask { CommonParameters = new CommonParameters(digestionParams: rnaDigestion) }, "Search");

        settings.Modifications.ResetToDefaults(settings.DigestionParametersForDefaults);

        var rnaDefaults = new CommonParameters(digestionParams: rnaDigestion);
        Assert.That(settings.Modifications.VariableSelection, Is.EquivalentTo(rnaDefaults.ListOfModsVariable));
        Assert.That(settings.Modifications.VariableSelection.Select(m => m.Item2),
            Does.Contain("Cyclic Phosphate on X"), "the RNA default, not Oxidation on M");
    }
}

/// <summary>
/// GPTMD's ListOfModsGptmd is the one setting users routinely change on a GPTMD task, and it uses the
/// same widget as the fixed and variable pickers.
/// </summary>
public class GptmdModificationTests
{
    [Test]
    public void ADefaultGptmdTaskLoadsItsModificationsAsSelected()
    {
        var settings = new TaskSettingsViewModel(new GptmdTask(), "GPTMD");

        Assert.That(settings.IsGptmdTask, Is.True);
        Assert.That(settings.Modifications.SelectedGptmdCount, Is.GreaterThan(0),
            "GptmdParameters ships a non-empty default list, which should arrive ticked");
    }

    [Test]
    public void GptmdModificationsRoundTrip()
    {
        var task = new GptmdTask();
        var settings = new TaskSettingsViewModel(task, "GPTMD");

        ModificationChoice spare = settings.Modifications.Groups.SelectMany(g => g.Choices)
            .First(c => !c.IsGptmd);
        spare.IsGptmd = true;
        settings.Apply();

        Assert.That(task.GptmdParameters.ListOfModsGptmd.Select(m => m.Item2), Does.Contain(spare.Name));

        spare.IsGptmd = false;
        settings.Apply();

        Assert.That(task.GptmdParameters.ListOfModsGptmd.Select(m => m.Item2), Does.Not.Contain(spare.Name),
            "unticking must remove it as well");
    }

    /// <summary>GPTMD is not exclusive with fixed or variable - it asks what to search for.</summary>
    [Test]
    public void GptmdSelectionIsIndependentOfFixedAndVariable()
    {
        var settings = new TaskSettingsViewModel(new GptmdTask(), "GPTMD");
        ModificationChoice choice = settings.Modifications.Groups.SelectMany(g => g.Choices).First();

        choice.IsGptmd = true;
        choice.IsVariable = true;

        Assert.That(choice.IsGptmd, Is.True, "choosing variable must not clear the GPTMD tick");
    }
}

/// <summary>The three task types added after Search, Calibrate and GPTMD.</summary>
public class RemainingTaskTypeTests
{
    [Test]
    public void CrosslinkSearchOptionsRoundTrip()
    {
        var task = new XLSearchTask();
        var settings = new TaskSettingsViewModel(task, "XLSearch");
        Assert.That(settings.IsXlSearchTask, Is.True);
        Assert.That(settings.IsSearchTask, Is.False, "an XL search is not a plain search");

        string otherCrosslinker = settings.Crosslinkers.First(c => c != settings.Crosslinker);
        settings.Crosslinker = otherCrosslinker;
        settings.CrosslinkSearchTopNum = 111;
        settings.CrosslinkAtCleavageSite = true;
        settings.Apply();

        Assert.That(task.XlSearchParameters.CrosslinkSearchTopNum, Is.EqualTo(111));
        Assert.That(task.XlSearchParameters.CrosslinkAtCleavageSite, Is.True);
        Assert.That(task.XlSearchParameters.Crosslinker.CrosslinkerName, Is.EqualTo(otherCrosslinker),
            "the crosslinker is chosen by name and looked up, so the object must actually change");
    }

    [Test]
    public void GlycoSearchOptionsRoundTrip()
    {
        var task = new GlycoSearchTask();
        var settings = new TaskSettingsViewModel(task, "GlycoSearch");
        Assert.That(settings.IsGlycoSearchTask, Is.True);

        settings.GlycoSearchType = GlycoSearchType.NGlycanSearch;
        settings.GlycoSearchTopNum = 27;
        settings.MaximumOGlycanAllowed = 6;
        settings.OxoniumIonFilt = false;
        settings.Apply();

        Assert.That(task._glycoSearchParameters.GlycoSearchType, Is.EqualTo(GlycoSearchType.NGlycanSearch));
        Assert.That(task._glycoSearchParameters.GlycoSearchTopNum, Is.EqualTo(27));
        Assert.That(task._glycoSearchParameters.MaximumOGlycanAllowed, Is.EqualTo(6));
        Assert.That(task._glycoSearchParameters.OxoniumIonFilt, Is.False);
    }

    [Test]
    public void SpectralAveragingOptionsRoundTrip()
    {
        var task = new SpectralAveragingTask();
        var settings = new TaskSettingsViewModel(task, "Average");
        Assert.That(settings.IsAveragingTask, Is.True);

        settings.NumberOfScansToAverage = 9;
        settings.ScanOverlap = 3;
        settings.Apply();

        Assert.That(task.Parameters.NumberOfScansToAverage, Is.EqualTo(9));
        Assert.That(task.Parameters.ScanOverlap, Is.EqualTo(3));
    }

    /// <summary>All six task types share the common settings, so all six must accept them.</summary>
    [Test]
    public void EverySupportedTaskTypeAcceptsCommonSettings()
    {
        MetaMorpheusTask[] tasks =
        {
            new SearchTask(), new CalibrationTask(), new GptmdTask(),
            new XLSearchTask(), new GlycoSearchTask(), new SpectralAveragingTask(),
        };

        foreach (MetaMorpheusTask task in tasks)
        {
            var settings = new TaskSettingsViewModel(task, task.GetType().Name)
            {
                PrecursorTolerance = "7",
                MaxMissedCleavages = 3,
            };
            settings.Apply();

            Assert.That(task.CommonParameters.PrecursorMassTolerance.Value, Is.EqualTo(7).Within(1e-9),
                $"{task.GetType().Name} should accept the precursor tolerance");
            Assert.That(task.CommonParameters.DigestionParams.MaxMissedCleavages, Is.EqualTo(3),
                $"{task.GetType().Name} should accept the digestion settings");
        }
    }
}
