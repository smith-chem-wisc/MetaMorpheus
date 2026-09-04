using System;
using Avalonia;
using Avalonia.Controls.ApplicationLifetimes;
using Avalonia.Markup.Xaml;
using MetaMorpheusAvalonia.ViewModels;
using MetaMorpheusAvalonia.Views;
using System.Diagnostics.CodeAnalysis;

namespace MetaMorpheusAvalonia;

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
public partial class App : Application
{
    public override void Initialize() => AvaloniaXamlLoader.Load(this);

    public override void OnFrameworkInitializationCompleted()
    {
        if (ApplicationLifetime is IClassicDesktopStyleApplicationLifetime desktop)
        {
            desktop.MainWindow = new MainWindow { DataContext = BuildViewModel() };
        }
        base.OnFrameworkInitializationCompleted();
    }

    /// <summary>
    /// The view model's constructor loads the modification, protease and glycan databases off disk, so
    /// a missing data directory or a malformed file throws here. Failing that way would kill the
    /// process before any window existed - no message, on the platforms least likely to have the
    /// layout right. Show the window with the reason in its log instead.
    /// </summary>
    private static MainWindowViewModel BuildViewModel()
    {
        try
        {
            return new MainWindowViewModel();
        }
        catch (Exception exception)
        {
            return MainWindowViewModel.ThatFailedToStart(exception);
        }
    }
}
