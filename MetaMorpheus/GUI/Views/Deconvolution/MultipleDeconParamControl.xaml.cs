using System.Windows;
using System.Windows.Controls;
using GuiFunctions;

namespace MetaMorpheusGUI;

public partial class MultipleDeconParamControl : UserControl
{
    public MultipleDeconParamControl()
    {
        InitializeComponent();
    }

    private void RemoveSubType_Click(object sender, RoutedEventArgs e)
    {
        if (sender is Button button && button.Tag is DeconParamsViewModel vm
            && DataContext is MultipleDeconParamsViewModel hostVm)
        {
            hostVm.RemoveSubType(vm);
        }
    }

    private void AddSubType_Click(object sender, RoutedEventArgs e)
    {
        if (DataContext is MultipleDeconParamsViewModel hostVm)
        {
            hostVm.AddSubType(hostVm.SelectedAddType);
        }
    }
}
