using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using MetaMorpheus.Avalonia.ViewModels;

namespace MetaMorpheus.Avalonia.Views;

/// <summary>Only the parts that genuinely need a window: file pickers and opening a browser.</summary>
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
                new FilePickerFileType("Spectra") { Patterns = new[] { "*.mzML", "*.mgf", "*.raw" } },
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
                new FilePickerFileType("Databases") { Patterns = new[] { "*.xml", "*.xml.gz", "*.fasta", "*.fa" } },
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
