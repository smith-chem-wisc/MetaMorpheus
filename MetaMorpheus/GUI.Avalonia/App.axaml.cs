using System;
using Avalonia;
using Avalonia.Controls.ApplicationLifetimes;
using Avalonia.Markup.Xaml;
using MetaMorpheus.Avalonia.ViewModels;
using MetaMorpheus.Avalonia.Views;

namespace MetaMorpheus.Avalonia;

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
