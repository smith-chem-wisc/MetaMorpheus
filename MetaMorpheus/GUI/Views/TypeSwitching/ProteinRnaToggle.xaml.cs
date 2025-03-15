using System.Windows.Controls;

namespace MetaMorpheusGUI;

/// <summary>
/// Interaction logic for ProteinRnaToggle.xaml
/// </summary>
public partial class ProteinRnaToggle : UserControl
{
    public ProteinRnaToggle()
    {
        InitializeComponent();
    }

    private void ToggleCheckBox_Checked(object sender, System.Windows.RoutedEventArgs e)
    {
        var checkBox = (CheckBox)sender;
        var checkState = checkBox.IsChecked;
        UpdateGUISettings.Globals.IsRnaMode = checkState ?? true;
        UpdateGUISettings.NotifyGlobalsChanged();
    }
}