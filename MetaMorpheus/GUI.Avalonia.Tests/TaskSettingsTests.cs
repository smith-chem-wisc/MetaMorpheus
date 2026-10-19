using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using Avalonia.Headless.NUnit;
using EngineLayer;
using MassSpectrometry;
using MetaMorpheusAvalonia.ViewModels;
using MetaMorpheusAvalonia.Views;
using MzLibUtil;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using TaskLayer;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;

namespace Test.AvaloniaGui;

/// <summary>
/// The settings dialog is only useful if what it collects reaches the task, and only safe if what it
/// does not collect survives. CommonParameters has private setters, so Apply() clones and overrides
/// rather than calling the constructor - which would default the settings this dialog never shows.
/// These assert the round trip rather than the rendering.
/// </summary>
public class TaskSettingsTests
{
    /// <summary>
    /// The dialog shows 13 of the settings CommonParameters holds. The other ~30 have to come through
    /// Apply() untouched, which the constructor could not do: passing a subset of its arguments leaves
    /// the rest at their defaults, so editing a tolerance would quietly reset deconvolution, windowing
    /// and separation type. Several of these change search results, so the failure would be invisible
    /// in the GUI and visible only in the output.
    /// </summary>
    [Test]
    public void ApplyPreservesSettingsTheDialogDoesNotShow()
    {
        var task = new SearchTask();

        // values differing from the constructor defaults, on settings this view model never shows
        task.CommonParameters = new CommonParameters(
            totalPartitions: 7,
            numberOfPeaksToKeepPerWindow: 42,
            addCompIons: true,
            assumeOrphanPeaksAreZ1Fragments: false,
            separationType: "RPLC");

        var settings = new TaskSettingsViewModel(task, "Search") { PrecursorTolerance = "17" };
        settings.Apply();

        Assert.That(task.CommonParameters.PrecursorMassTolerance.Value, Is.EqualTo(17).Within(1e-9),
            "the edited setting still has to arrive");

        Assert.Multiple(() =>
        {
            Assert.That(task.CommonParameters.TotalPartitions, Is.EqualTo(7),
                "TotalPartitions was reset by Apply()");
            Assert.That(task.CommonParameters.NumberOfPeaksToKeepPerWindow, Is.EqualTo(42),
                "NumberOfPeaksToKeepPerWindow was reset by Apply()");
            Assert.That(task.CommonParameters.AddCompIons, Is.True,
                "AddCompIons was reset by Apply()");
            Assert.That(task.CommonParameters.AssumeOrphanPeaksAreZ1Fragments, Is.False,
                "AssumeOrphanPeaksAreZ1Fragments was reset by Apply()");
            Assert.That(task.CommonParameters.SeparationType, Is.EqualTo("RPLC"),
                "SeparationType was reset by Apply()");
        });
    }

    [Test]
    public void SettingsLoadFromTheTaskTheyWillEdit()
    {
        var task = new SearchTask();
        var settings = new TaskSettingsViewModel(task, "Search");

        Assert.That(settings.IsSearchTask, Is.True);
        Assert.That(settings.PrecursorTolerance, Is.EqualTo(
            ((PpmTolerance)task.CommonParameters.PrecursorMassTolerance).Value
                .ToString(System.Globalization.CultureInfo.InvariantCulture)));
        Assert.That(settings.DecoyType, Is.EqualTo(task.SearchParameters.DecoyType));
    }

    /// <summary>
    /// Tolerances are the setting most likely to be silently lost, since Apply() rebuilds
    /// CommonParameters and a missed named argument would quietly restore the default of 5 ppm.
    /// </summary>
    [Test]
    public void EditedTolerancesReachTheTask()
    {
        var task = new SearchTask();
        var settings = new TaskSettingsViewModel(task, "Search")
        {
            PrecursorTolerance = "13",
            PrecursorToleranceIsPpm = true,
            ProductTolerance = "0.4",
            ProductToleranceIsPpm = false,
        };

        settings.Apply();

        Assert.That(task.CommonParameters.PrecursorMassTolerance, Is.TypeOf<PpmTolerance>());
        Assert.That(task.CommonParameters.PrecursorMassTolerance.Value, Is.EqualTo(13).Within(1e-9));
        Assert.That(task.CommonParameters.ProductMassTolerance, Is.TypeOf<AbsoluteTolerance>(),
            "switching to absolute must change the tolerance type, not only the number");
        Assert.That(task.CommonParameters.ProductMassTolerance.Value, Is.EqualTo(0.4).Within(1e-9));
    }

    [Test]
    public void EditedDigestionSettingsReachTheTask()
    {
        var task = new SearchTask();
        var settings = new TaskSettingsViewModel(task, "Search")
        {
            MaxMissedCleavages = 4,
            MinPeptideLength = 6,
            MaxPeptideLength = 40,
            MaxModsPerPeptide = 3,
        };

        settings.Apply();

        Assert.That(task.CommonParameters.DigestionParams.MaxMissedCleavages, Is.EqualTo(4));
        Assert.That(task.CommonParameters.DigestionParams.MinLength, Is.EqualTo(6));
        Assert.That(task.CommonParameters.DigestionParams.MaxLength, Is.EqualTo(40));
        Assert.That(task.CommonParameters.DigestionParams.MaxMods, Is.EqualTo(3));
    }

    [Test]
    public void EditedSearchOptionsReachTheTask()
    {
        var task = new SearchTask();
        var settings = new TaskSettingsViewModel(task, "Search")
        {
            DoParsimony = false,
            NoOneHitWonders = true,
            MatchBetweenRuns = true,
            WritePrunedDatabase = true,
            DecoyType = DecoyType.None,
            MassDiffAcceptorType = MassDiffAcceptorType.Exact,
        };

        settings.Apply();

        Assert.That(task.SearchParameters.DoParsimony, Is.False);
        Assert.That(task.SearchParameters.NoOneHitWonders, Is.True);
        Assert.That(task.SearchParameters.MatchBetweenRuns, Is.True);
        Assert.That(task.SearchParameters.WritePrunedDatabase, Is.True);
        Assert.That(task.SearchParameters.DecoyType, Is.EqualTo(DecoyType.None));
        Assert.That(task.SearchParameters.MassDiffAcceptorType, Is.EqualTo(MassDiffAcceptorType.Exact));
    }

    /// <summary>
    /// Apply() rebuilds CommonParameters from a subset of its 42 constructor arguments, so anything
    /// not represented in the dialog must survive rather than reverting to a constructor default.
    /// </summary>
    [Test]
    public void SettingsNotShownInTheDialogAreNotLost()
    {
        var task = new SearchTask
        {
            CommonParameters = new CommonParameters(
                taskDescriptor: "my-descriptor",
                listOfModsFixed: new[] { ("Common Fixed", "Carbamidomethyl on C") }),
        };

        new TaskSettingsViewModel(task, "Search").Apply();

        Assert.That(task.CommonParameters.TaskDescriptor, Is.EqualTo("my-descriptor"),
            "the task descriptor is not editable here, so it must be carried through Apply");
        Assert.That(task.CommonParameters.ListOfModsFixed.Select(m => m.Item2),
            Does.Contain("Carbamidomethyl on C"),
            "fixed modifications are not editable here yet, so they must survive Apply");
    }

    [Test]
    public void DissociationTypeReachesTheTask()
    {
        var task = new SearchTask();
        var settings = new TaskSettingsViewModel(task, "Search") { DissociationType = DissociationType.ETD };

        settings.Apply();

        Assert.That(task.CommonParameters.DissociationType, Is.EqualTo(DissociationType.ETD));
    }

    [TestCase("0", "precursor")]
    [TestCase("-1", "precursor")]
    [TestCase("abc", "precursor")]
    public void NonsenseTolerancesAreReportedRatherThanApplied(string value, string expectedMention)
    {
        var settings = new TaskSettingsViewModel(new SearchTask(), "Search") { PrecursorTolerance = value };

        var problems = settings.Validate();

        Assert.That(problems, Is.Not.Empty);
        Assert.That(string.Join(" ", problems).ToLowerInvariant(), Does.Contain(expectedMention));
    }

    [Test]
    public void PeptideLengthRangeIsValidated()
    {
        var settings = new TaskSettingsViewModel(new SearchTask(), "Search")
        {
            MinPeptideLength = 20,
            MaxPeptideLength = 5,
        };

        Assert.That(settings.Validate(), Is.Not.Empty);
    }

    [Test]
    public void ValidSettingsProduceNoComplaints()
    {
        Assert.That(new TaskSettingsViewModel(new SearchTask(), "Search").Validate(), Is.Empty);
    }

    /// <summary>Non-search tasks share the common settings and hide the search-only section.</summary>
    [Test]
    public void CalibrationAndGptmdTasksAreEditableToo()
    {
        foreach (MetaMorpheusTask task in new MetaMorpheusTask[] { new CalibrationTask(), new GptmdTask() })
        {
            var settings = new TaskSettingsViewModel(task, "x") { PrecursorTolerance = "9" };
            settings.Apply();

            Assert.That(settings.IsSearchTask, Is.False, $"{task.GetType().Name} is not a search");
            Assert.That(task.CommonParameters.PrecursorMassTolerance.Value, Is.EqualTo(9).Within(1e-9),
                $"{task.GetType().Name} should accept common settings");
        }
    }

    /// <summary>
    /// Gate A. Load a task holding non-default values across the parameter space, click Save without
    /// editing anything, and every setting must come back unchanged. Reflective over CommonParameters,
    /// so a property added later is covered without this test having to know about it.
    /// </summary>
    [Test]
    public void ApplyWithoutEditingChangesNothing()
    {
        var task = new SearchTask
        {
            CommonParameters = new CommonParameters(
                taskDescriptor: "descriptor",
                ms2childScanDissociationType: DissociationType.EThcD,
                separationType: "CZE",
                addCompIons: true,
                totalPartitions: 3,
                numberOfPeaksToKeepPerWindow: 1000,
                trimMs1Peaks: true,
                addTruncations: true,
                pepQValueThreshold: 0.02,
                precursorMassMatchMode: PrecursorMassMatchMode.MostAbundant,
                digestionParams: new DigestionParams(
                    protease: "Asp-N",
                    initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                    searchModeType: CleavageSpecificity.Semi,
                    fragmentationTerminus: FragmentationTerminus.C,
                    keepNGlycopeptide: true,
                    keepOGlycopeptide: true)),
        };

        Dictionary<string, string> before = Snapshot(task.CommonParameters);

        new TaskSettingsViewModel(task, "Search").Apply();

        Assert.That(Snapshot(task.CommonParameters), Is.EqualTo(before));
    }

    /// <summary>The digestion half of the same contract, spelled out so a failure names the setting.</summary>
    [Test]
    public void DigestionSettingsNotShownInTheDialogAreNotLost()
    {
        var task = new SearchTask
        {
            CommonParameters = new CommonParameters(digestionParams: new DigestionParams(
                protease: "Asp-N",
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                searchModeType: CleavageSpecificity.Semi,
                fragmentationTerminus: FragmentationTerminus.C,
                generateUnlabeledProteinsForSilac: false,
                keepNGlycopeptide: true,
                keepOGlycopeptide: true)),
        };

        // edit something the dialog does show, so this is a realistic Save rather than a no-op
        new TaskSettingsViewModel(task, "Search") { MaxMissedCleavages = 4 }.Apply();

        var digestion = (DigestionParams)task.CommonParameters.DigestionParams;
        Assert.Multiple(() =>
        {
            Assert.That(digestion.MaxMissedCleavages, Is.EqualTo(4), "the edited value should arrive");
            Assert.That(digestion.SearchModeType, Is.EqualTo(CleavageSpecificity.Semi));
            Assert.That(digestion.FragmentationTerminus, Is.EqualTo(FragmentationTerminus.C));
            Assert.That(digestion.InitiatorMethionineBehavior, Is.EqualTo(InitiatorMethionineBehavior.Retain));
            Assert.That(digestion.GeneratehUnlabeledProteinsForSilac, Is.False);
            Assert.That(digestion.KeepNGlycopeptide, Is.True);
            Assert.That(digestion.KeepOGlycopeptide, Is.True);
        });
    }

    /// <summary>
    /// A non-specific search keeps the user's enzyme in SpecificProtease and singleN/singleC in
    /// Protease, so the dialog has to read and write the former.
    /// </summary>
    [Test]
    public void NonSpecificSearchKeepsItsModeAndShowsTheSpecificProtease()
    {
        var task = new SearchTask
        {
            CommonParameters = new CommonParameters(digestionParams: new DigestionParams(
                protease: "Asp-N",
                searchModeType: CleavageSpecificity.None,
                fragmentationTerminus: FragmentationTerminus.N)),
        };

        var settings = new TaskSettingsViewModel(task, "Search");
        Assert.That(settings.Protease, Is.EqualTo("Asp-N"), "the enzyme box should show the user's choice, not singleN");

        settings.Apply();

        var digestion = (DigestionParams)task.CommonParameters.DigestionParams;
        Assert.Multiple(() =>
        {
            Assert.That(digestion.SearchModeType, Is.EqualTo(CleavageSpecificity.None), "the search mode must survive Save");
            Assert.That(digestion.SpecificProtease.Name, Is.EqualTo("Asp-N"));
            Assert.That(digestion.Protease.Name, Is.EqualTo("singleN"));
        });
    }

    /// <summary>An RNA task must stay an RNA task; the enzyme box has nothing to offer it.</summary>
    [Test]
    public void RnaDigestionParametersSurviveApply()
    {
        var task = new SearchTask
        {
            CommonParameters = new CommonParameters(
                digestionParams: new RnaDigestionParams(rnase: "RNase T1", maxMissedCleavages: 1)),
        };

        new TaskSettingsViewModel(task, "Search") { MinPeptideLength = 4 }.Apply();

        Assert.That(task.CommonParameters.DigestionParams, Is.TypeOf<RnaDigestionParams>());
        var rna = (RnaDigestionParams)task.CommonParameters.DigestionParams;
        Assert.Multiple(() =>
        {
            Assert.That(rna.Rnase.Name, Is.EqualTo("RNase T1"));
            Assert.That(rna.MaxMissedCleavages, Is.EqualTo(1));
            Assert.That(rna.MinLength, Is.EqualTo(4), "the edited value should still arrive");
        });
    }

    /// <summary>
    /// On de-DE, "." is the group separator, so a current-culture double.TryParse("0.05") returns true
    /// with the value 5.0 - a tolerance a hundred times looser than the one typed, with no error shown.
    /// This is the locale bug the PR most needed not to have, given what it is for.
    /// </summary>
    [TestCase("de-DE")]
    [TestCase("fr-FR")]
    [TestCase("en-US")]
    public void TolerancesRoundTripUnderAnyCulture(string culture)
    {
        CultureInfo original = CultureInfo.CurrentCulture;
        try
        {
            CultureInfo.CurrentCulture = new CultureInfo(culture);

            var task = new SearchTask();
            var settings = new TaskSettingsViewModel(task, "Search")
            {
                ProductTolerance = "0.05",
                ProductToleranceIsPpm = false,
            };
            settings.Apply();

            Assert.That(task.CommonParameters.ProductMassTolerance.Value, Is.EqualTo(0.05).Within(1e-12),
                $"0.05 was misread under {culture}");

            // and back out again, so a value written on one locale is not misread on another
            Assert.That(new TaskSettingsViewModel(task, "Search").ProductTolerance, Is.EqualTo("0.05"));
        }
        finally
        {
            CultureInfo.CurrentCulture = original;
        }
    }

    /// <summary>Same field, so the two GUIs must not disagree about what is valid.</summary>
    [TestCase(0.5, false)]
    [TestCase(1, true)]
    [TestCase(5, true)]
    public void ScoreCutoffMatchesTaskValidatorsThreshold(double cutoff, bool valid)
    {
        var settings = new TaskSettingsViewModel(new SearchTask(), "Search") { ScoreCutoff = cutoff };

        Assert.That(settings.Validate(), valid ? Is.Empty : Is.Not.Empty);
    }

    private static Dictionary<string, string> Snapshot(CommonParameters parameters) =>
        typeof(CommonParameters).GetProperties()
            .Where(p => p.CanRead && p.GetIndexParameters().Length == 0)
            .ToDictionary(p => p.Name, p => Describe(p.GetValue(parameters)));

    private static string Describe(object value) => value switch
    {
        null => "<null>",
        IEnumerable<(string, string)> mods => string.Join(";", mods.Select(m => $"{m.Item1}|{m.Item2}")),
        _ => value.ToString(),
    };

    [AvaloniaTest]
    public void SettingsWindowBuilds()
    {
        var window = new TaskSettingsWindow
        {
            DataContext = new TaskSettingsViewModel(new SearchTask(), "Search"),
        };
        window.Show();

        Assert.That(window.IsVisible, Is.True);
        Assert.That(window.Saved, Is.False, "nothing is saved until the user clicks Save");
    }
}
