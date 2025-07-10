using GuiFunctions;
using System;
using System.Windows;
using System.Windows.Controls;

namespace MetaMorpheusGUI;

/// <summary>
/// Interaction logic for SettingsButtonControl.xaml
/// </summary>
public partial class SettingsButtonControl : UserControl
{
    public SettingsButtonControl()
    {
        InitializeComponent();
        var vm = MetaDrawSettingsViewModel.Instance; 
        SettingsViewModel = vm; // Assign the view model to the control's property
    }

    // DependencyProperty for the settings view model
    public static readonly DependencyProperty SettingsViewModelProperty =
        DependencyProperty.Register(nameof(SettingsViewModel), typeof(MetaDrawSettingsViewModel), typeof(SettingsButtonControl));

    public MetaDrawSettingsViewModel SettingsViewModel
    {
        get => GetValue(SettingsViewModelProperty) as MetaDrawSettingsViewModel;
        set => SetValue(SettingsViewModelProperty, value);
    }

    // DependencyProperty for the refresh action
    public static readonly DependencyProperty RefreshActionProperty =
        DependencyProperty.Register(nameof(RefreshAction), typeof(Action), typeof(SettingsButtonControl));

    public Action RefreshAction
    {
        get => GetValue(RefreshActionProperty) as Action;
        set => SetValue(RefreshActionProperty, value);
    }

    private void SettingsButton_Click(object sender, RoutedEventArgs e)
    {
        // Open the settings window
        var settingsWindow = new MetaDrawSettingsWindow()
        {
            Owner = Window.GetWindow(this)
        };

        var result = settingsWindow.ShowDialog();

        // If settings were changed, call the refresh action
        if (result == true)
        {
            RefreshAction?.Invoke();
        }
    }
}
