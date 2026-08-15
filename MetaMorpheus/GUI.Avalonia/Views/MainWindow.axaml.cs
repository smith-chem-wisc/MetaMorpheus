using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using EngineLayer;
using MetaMorpheus.Avalonia.ViewModels;

namespace MetaMorpheus.Avalonia.Views;

/// <summary>Only the parts that genuinely need a window: file pickers and opening a browser.</summary>
public partial class MainWindow : Window
{
    private MainWindowViewModel ViewModel => (MainWindowViewModel)DataContext;

    public MainWindow() => AvaloniaXamlLoader.Load(this);

    /// <summary>
    /// Built from GlobalVariables rather than written out, so the picker cannot drift from what the
    /// application accepts. Two things this gets right that a hand-written list did not:
    ///
    ///   * The entries are lowercase. Avalonia's FreeDesktop backend passes these to the XDG portal as
    ///     GlobStyle globs, and portal matching is case-sensitive, so "*.mzML" hides sample.mzml on
    ///     Linux - the platform this exists for. Windows and macOS match case-insensitively, so it
    ///     would never reproduce locally on either.
    ///   * ".d" is excluded. Bruker data is a directory and OpenFilePickerAsync cannot select one;
    ///     users pick the .tdf inside it and AddSpectraFiles maps that back to the folder.
    /// </summary>
    private static string[] SpectraPatterns => GlobalVariables.AcceptedSpectraFormats
        .Where(extension => extension != ".d")
        .Select(extension => "*" + extension)
        .ToArray();

    private static string[] DatabasePatterns => GlobalVariables.AcceptedDatabaseFormats
        .SelectMany(extension => new[] { "*" + extension, "*" + extension + ".gz" })
        .ToArray();

    private async void OnAddSpectraClick(object sender, RoutedEventArgs e)
    {
        IReadOnlyList<IStorageFile> files = await StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Select spectra files",
            AllowMultiple = true,
            FileTypeFilter = new[]
            {
                new FilePickerFileType("Spectra") { Patterns = SpectraPatterns },
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
                new FilePickerFileType("Databases") { Patterns = DatabasePatterns },
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
