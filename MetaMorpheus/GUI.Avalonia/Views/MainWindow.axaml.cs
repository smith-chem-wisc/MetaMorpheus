using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Diagnostics;
using System.Linq;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using EngineLayer;
using MetaMorpheusAvalonia.ViewModels;

namespace MetaMorpheusAvalonia.Views;

/// <summary>Only the parts that genuinely need a window: file pickers and opening a browser.</summary>
// Excluded from coverage deliberately, matching how the WPF GUI is treated.
//
// PRECEDENT: Test.csproj references CMD, EngineLayer, GuiFunctions and TaskLayer - but NOT
// GUI.csproj. The WPF views are therefore never loaded by the test host, never instrumented, and
// contribute nothing to the coverage denominator: Codecov reports 0 files under MetaMorpheus/GUI/
// and 79 under MetaMorpheus/GuiFunctions/. That split is the whole point of GuiFunctions - the
// logic was moved out of the views so that it COULD be covered, leaving behind only the parts that
// cannot meaningfully be unit tested.
//
// GUI.Avalonia cannot use the same mechanism, because its tests construct real windows headlessly
// and so must reference this project. [ExcludeFromCodeCoverage] is the equivalent, and is already
// the house instrument for untestable code (~20 files on master carry it, e.g. the design-time
// view models in GuiFunctions/ViewModels/Deconvolution).
//
// THE RULE THIS ENCODES: if something here can be tested, it does not belong here - it belongs on
// the view model, where it is covered like the rest of GuiFunctions. That is why SpectraPatterns
// and DatabasePatterns now live on MainWindowViewModel rather than in this file. Do not add logic
// to a class carrying this attribute in order to avoid writing a test for it.
[ExcludeFromCodeCoverage]
public partial class MainWindow : Window
{
    private MainWindowViewModel ViewModel => (MainWindowViewModel)DataContext;

    public MainWindow() => AvaloniaXamlLoader.Load(this);

    private async void OnAddSpectraClick(object sender, RoutedEventArgs e)
    {
        IReadOnlyList<IStorageFile> files = await StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Select spectra files",
            AllowMultiple = true,
            FileTypeFilter = new[]
            {
                new FilePickerFileType("Spectra") { Patterns = MainWindowViewModel.SpectraPatterns },
            },
        });
        ViewModel.AddSpectraFiles(files.Select(f => f.Path.LocalPath));
    }

    private async void OnAddDatabaseClick(object sender, RoutedEventArgs e)
    {
        IReadOnlyList<IStorageFile> files = await StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Select protein databases",
            AllowMultiple = true,
            FileTypeFilter = new[]
            {
                new FilePickerFileType("Databases") { Patterns = MainWindowViewModel.DatabasePatterns },
            },
        });
        ViewModel.AddDatabases(files.Select(f => f.Path.LocalPath));
    }

    private async void OnBrowseOutputClick(object sender, RoutedEventArgs e)
    {
        IReadOnlyList<IStorageFolder> folders = await StorageProvider.OpenFolderPickerAsync(
            new FolderPickerOpenOptions { Title = "Select the output folder", AllowMultiple = false });
        if (folders.Count > 0)
        {
            ViewModel.OutputFolder = folders[0].Path.LocalPath;
        }
    }

    /// <summary>Opens the settings dialog for the selected task and refreshes the list if it saved.</summary>
    private async void OnEditTaskClick(object sender, RoutedEventArgs e)
    {
        var selected = this.FindControl<ListBox>("TaskList")?.SelectedItem;
        if (selected is null)
        {
            ViewModel.Note("Select a task first, then edit its settings.");
            return;
        }

        var settings = ViewModel.CreateSettingsFor(selected);
        if (settings is null)
        {
            // Reachable in a way the WPF equivalent is not: this comes from a ListBox item rather than
            // an exhaustive enum, so an unexpected item type must not take the window down.
            ViewModel.Note("That item is not an editable task.");
            return;
        }

        var dialog = new TaskSettingsWindow { DataContext = settings };
        await dialog.ShowDialog(this);
        if (dialog.Saved)
        {
            ViewModel.Note($"Updated settings for {settings.TaskKind}.");
        }
    }

    private void OnWikiClick(object sender, RoutedEventArgs e) =>
        OpenUrl("https://github.com/smith-chem-wisc/MetaMorpheus/wiki");

    private void OnIssuesClick(object sender, RoutedEventArgs e) =>
        OpenUrl("https://github.com/smith-chem-wisc/MetaMorpheus/issues");

    private static void OpenUrl(string url)
    {
        try { Process.Start(new ProcessStartInfo(url) { UseShellExecute = true }); }
        catch { /* opening a browser is a convenience, not worth taking the window down for */ }
    }
}
