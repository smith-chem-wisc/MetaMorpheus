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
    public event EventHandler<MetaDrawSettingsChangedEventArgs> SettingsChanged;

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

    public static readonly DependencyProperty SelectedTabIndexProperty =
        DependencyProperty.Register(
            nameof(SelectedTabIndex),
            typeof(int),
            typeof(MetaDrawSettingsWindow),
            new PropertyMetadata(0));

    public int SelectedTabIndex
    {
        get => (int)GetValue(SelectedTabIndexProperty);
        set => SetValue(SelectedTabIndexProperty, value);
    }

    private void SettingsButton_Click(object sender, RoutedEventArgs e)
    {
        // Open the settings window
        var settingsWindow = new MetaDrawSettingsWindow()
        {
            Owner = Window.GetWindow(this)
        };

        settingsWindow.SettingsWindowTabControl.SelectedIndex = SelectedTabIndex;

        var originalQ = MetaDrawSettings.QValueFilter;
        var originalShowDecoys = MetaDrawSettings.ShowDecoys;
        var originalAmbiguityFilter = MetaDrawSettings.AmbiguityFilter;
        var originalGlycoLevelMin = MetaDrawSettings.LocalizationLevelStart;
        var originalGlycoLevelMax = MetaDrawSettings.LocalizationLevelEnd;

        var result = settingsWindow.ShowDialog();

        // If canceled or closed, return and don't save settings. 
        if (result != true) 
            return;

        var args = new MetaDrawSettingsChangedEventArgs();

        // If any filtering setting has changed, set the filter changed flag to ensure the data grids are refreshed. 
        if (Math.Abs(MetaDrawSettings.QValueFilter - originalQ) > 0.00001 
            || originalShowDecoys != MetaDrawSettings.ShowDecoys
            || originalAmbiguityFilter != MetaDrawSettings.AmbiguityFilter
            || originalGlycoLevelMin != MetaDrawSettings.LocalizationLevelStart
            || originalGlycoLevelMax != MetaDrawSettings.LocalizationLevelEnd
            )
        {
            args.FilterChanged = true;
        }
        SettingsChanged?.Invoke(this, args);
    }
}

public class MetaDrawSettingsChangedEventArgs : EventArgs
{
    public bool FilterChanged { get; set; } = false;
}
