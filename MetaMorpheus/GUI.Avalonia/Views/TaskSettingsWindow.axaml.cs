using System.Linq;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using MetaMorpheusAvalonia.ViewModels;

namespace MetaMorpheusAvalonia.Views;

public partial class TaskSettingsWindow : Window
{
    private TaskSettingsViewModel ViewModel => (TaskSettingsViewModel)DataContext;

    public TaskSettingsWindow() => AvaloniaXamlLoader.Load(this);

    /// <summary>True when the user saved. The caller uses this to know whether to refresh.</summary>
    public bool Saved { get; private set; }

    private void OnSaveClick(object sender, RoutedEventArgs e)
    {
        var problems = ViewModel.Validate();
        if (problems.Count > 0)
        {
            // Report in place rather than in a second dialog, and do not close.
            this.FindControl<TextBlock>("ValidationText").Text = string.Join(" ", problems);
            return;
        }
        ViewModel.Apply();
        Saved = true;
        Close();
    }

    private void OnCancelClick(object sender, RoutedEventArgs e) => Close();

    /// <summary>Defaults come from CommonParameters, so RNA tasks get RNA's rather than the peptide set.</summary>
    private void OnResetModificationsClick(object sender, RoutedEventArgs e) =>
        ViewModel.Modifications.ResetToDefaults(ViewModel.DigestionParametersForDefaults);
}
