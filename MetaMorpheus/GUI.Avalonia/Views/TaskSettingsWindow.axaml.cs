using System.Linq;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using MetaMorpheusAvalonia.ViewModels;
using System.Diagnostics.CodeAnalysis;

namespace MetaMorpheusAvalonia.Views;

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
