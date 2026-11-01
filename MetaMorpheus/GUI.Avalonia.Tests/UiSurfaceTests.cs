using System.Linq;
using Avalonia.Controls;
using Avalonia.Headless.NUnit;
using Avalonia.VisualTree;
using EngineLayer;
using MetaMorpheusAvalonia.ViewModels;
using MetaMorpheusAvalonia.Views;
using NUnit.Framework;
using TaskLayer;

namespace Test.AvaloniaGui;

/// <summary>
/// Gate F: everything the PR description claims has to be reachable by clicking. These assert against
/// the real windows, because a capability that exists only on a view model is not a feature.
/// </summary>
public class UiSurfaceTests
{
    /// <summary>All six task types the description lists must actually be addable.</summary>
    [AvaloniaTest]
    public void AllSixTaskTypesCanBeAdded()
    {
        using var viewModel = new MainWindowViewModel();

        viewModel.AddCalibrationTaskCommand.Execute(null);
        viewModel.AddGptmdTaskCommand.Execute(null);
        viewModel.AddSearchTaskCommand.Execute(null);
        viewModel.AddXlSearchTaskCommand.Execute(null);
        viewModel.AddGlycoSearchTaskCommand.Execute(null);
        viewModel.AddAveragingTaskCommand.Execute(null);

        Assert.That(viewModel.Tasks.Select(t => t.Task.GetType()), Is.EquivalentTo(new[]
        {
            typeof(CalibrationTask), typeof(GptmdTask), typeof(SearchTask),
            typeof(XLSearchTask), typeof(GlycoSearchTask), typeof(SpectralAveragingTask),
        }));
    }

    /// <summary>
    /// Each task type's settings dialog must show its own section, or the Apply() arm behind it is
    /// unreachable code that only looks supported.
    /// </summary>
    [AvaloniaTest]
    [TestCase(typeof(SearchTask), nameof(TaskSettingsViewModel.IsSearchTask))]
    [TestCase(typeof(GptmdTask), nameof(TaskSettingsViewModel.IsGptmdTask))]
    [TestCase(typeof(XLSearchTask), nameof(TaskSettingsViewModel.IsXlSearchTask))]
    [TestCase(typeof(GlycoSearchTask), nameof(TaskSettingsViewModel.IsGlycoSearchTask))]
    [TestCase(typeof(SpectralAveragingTask), nameof(TaskSettingsViewModel.IsAveragingTask))]
    public void EachTaskTypesSettingsSectionIsShown(System.Type taskType, string flag)
    {
        var task = (MetaMorpheusTask)System.Activator.CreateInstance(taskType);
        var settings = new TaskSettingsViewModel(task, taskType.Name);

        var window = new TaskSettingsWindow { DataContext = settings };
        window.Show();

        bool shown = (bool)typeof(TaskSettingsViewModel).GetProperty(flag).GetValue(settings);
        Assert.That(shown, Is.True, $"{taskType.Name} should light up {flag}, which is what shows its section");
    }

    /// <summary>
    /// The modification picker has to exist in the dialog, not only on the view model. This also proves
    /// the compiled-binding source generator is running: x:Name resolution is what it produces, and it
    /// was being skipped while the SDK shipped an older Roslyn.
    /// </summary>
    [AvaloniaTest]
    public void TheModificationPickerIsInTheDialog()
    {
        var settings = new TaskSettingsViewModel(new SearchTask(), "Search");
        var window = new TaskSettingsWindow { DataContext = settings };
        window.Show();

        var list = window.FindControl<ItemsControl>("ModificationList");
        Assert.That(list, Is.Not.Null, "the modifications expander should contain the grouped list");
        Assert.That(list.ItemsSource, Is.SameAs(settings.Modifications.Groups));
        Assert.That(settings.Modifications.Groups, Is.Not.Empty);
    }

    [AvaloniaTest]
    public void TheGptmdModificationPickerIsInTheDialog()
    {
        var settings = new TaskSettingsViewModel(new GptmdTask(), "GPTMD");
        var window = new TaskSettingsWindow { DataContext = settings };
        window.Show();

        Assert.That(window.FindControl<ItemsControl>("GptmdModificationList"), Is.Not.Null);
        Assert.That(settings.IsGptmdTask, Is.True, "which is what makes the GPTMD section visible");
    }

    /// <summary>
    /// The contaminant flag is a filename guess, so the user must be able to overrule it and the run
    /// must honour the override. Previously the grid was read-only and whatever the heuristic decided
    /// was final.
    /// </summary>
    [Test]
    public void TheContaminantGuessCanBeOverridden()
    {
        using var viewModel = new MainWindowViewModel();
        viewModel.AddDatabases(new[] { "some-contaminants.fasta", "ordinary.fasta" });

        DatabaseForDisplay guessed = viewModel.Databases.Single(d => d.FileName.Contains("contaminants"));
        DatabaseForDisplay ordinary = viewModel.Databases.Single(d => d.FileName == "ordinary.fasta");

        Assert.Multiple(() =>
        {
            Assert.That(guessed.IsContaminant, Is.True, "the filename heuristic should have flagged it");
            Assert.That(ordinary.IsContaminant, Is.False);
        });

        guessed.IsContaminant = false;
        ordinary.IsContaminant = true;

        Assert.Multiple(() =>
        {
            Assert.That(guessed.IsContaminant, Is.False, "the user's correction must stick");
            Assert.That(ordinary.IsContaminant, Is.True);
        });
    }

    /// <summary>The grid column has to be writable, or the settable property is unreachable.</summary>
    [AvaloniaTest]
    public void TheContaminantColumnIsEditable()
    {
        using var viewModel = new MainWindowViewModel();
        var window = new MainWindow { DataContext = viewModel };
        window.Show();

        DataGrid grid = window.GetVisualDescendants().OfType<DataGrid>()
            .Single(g => g.Columns.Any(c => c.Header as string == "Contaminant"));

        Assert.Multiple(() =>
        {
            Assert.That(grid.IsReadOnly, Is.False, "a read-only grid makes every column read-only");
            Assert.That(grid.Columns.Single(c => c.Header as string == "Contaminant").IsReadOnly, Is.False);
            Assert.That(grid.Columns.Single(c => c.Header as string == "File").IsReadOnly, Is.True,
                "the filename is not something to edit");
        });
    }

    /// <summary>
    /// "+ADD DEFAULT CONTAMINANTS" exists in the WPF window and was missing here entirely. The
    /// contaminant list ships with the application, so this must find it wherever it is installed.
    /// </summary>
    [Test]
    public void DefaultContaminantsCanBeAdded()
    {
        using var viewModel = new MainWindowViewModel();
        viewModel.AddDefaultContaminantsCommand.Execute(null);

        Assert.That(viewModel.Databases, Is.Not.Empty,
            $"nothing loaded. Log:{System.Environment.NewLine}{viewModel.Log}");
        Assert.That(viewModel.Databases.All(d => d.IsContaminant), Is.True,
            "everything from the Contaminants folder should arrive already flagged");
    }

    /// <summary>
    /// Offering Custom in the combo box without carrying the expression meant the choice could be made
    /// and then silently ignored. The WPF window has always written it.
    /// </summary>
    [Test]
    public void ACustomMassDifferenceAcceptorReachesTheTask()
    {
        var task = new SearchTask();
        var settings = new TaskSettingsViewModel(task, "Search")
        {
            MassDiffAcceptorType = MassDiffAcceptorType.Custom,
            CustomMassDiffAcceptor = "myAcceptor dot 5 ppm 0,1.003,2.006",
        };

        Assert.That(settings.Validate(), Is.Empty, "that expression is valid, so it should not be rejected");
        settings.Apply();

        Assert.That(task.SearchParameters.MassDiffAcceptorType, Is.EqualTo(MassDiffAcceptorType.Custom));
        Assert.That(task.SearchParameters.CustomMdac, Is.EqualTo("myAcceptor dot 5 ppm 0,1.003,2.006"));

        // and it survives a reload, so the dialog shows what the task holds
        Assert.That(new TaskSettingsViewModel(task, "Search").CustomMassDiffAcceptor,
            Is.EqualTo("myAcceptor dot 5 ppm 0,1.003,2.006"));
    }

    [Test]
    public void AnUnreadableCustomAcceptorIsRejectedRatherThanRunning()
    {
        var settings = new TaskSettingsViewModel(new SearchTask(), "Search")
        {
            MassDiffAcceptorType = MassDiffAcceptorType.Custom,
            CustomMassDiffAcceptor = "this is not an acceptor",
        };

        Assert.That(settings.Validate(), Is.Not.Empty,
            "an unparseable expression should be caught here, not part-way into a search");
    }

    /// <summary>
    /// A ListBox hands back whatever object it holds, so an unexpected item type must not take the
    /// window down. The WPF equivalent switches on an exhaustive enum, so it cannot reach this state.
    /// </summary>
    [Test]
    public void AnItemThatIsNotATaskProducesNoSettings()
    {
        using var viewModel = new MainWindowViewModel();

        Assert.That(viewModel.CreateSettingsFor("not a task"), Is.Null);
        Assert.That(viewModel.CreateSettingsFor(null), Is.Null);
    }

    /// <summary>
    /// The log is appended to thousands of times during a search, so it is a StringBuilder rather than
    /// a string being concatenated. Appending from several threads at once must not lose lines or throw.
    /// </summary>
    [Test]
    public void TheLogSurvivesConcurrentAppends()
    {
        using var viewModel = new MainWindowViewModel();
        int before = viewModel.Log.Split(System.Environment.NewLine).Length;

        System.Threading.Tasks.Parallel.For(0, 200, i => viewModel.Note($"line-{i}"));

        string log = viewModel.Log;
        Assert.That(log.Split(System.Environment.NewLine).Length, Is.EqualTo(before + 200),
            "every line should be present exactly once");
        Assert.That(Enumerable.Range(0, 200).All(i => log.Contains($"line-{i}")), Is.True);
    }
}
