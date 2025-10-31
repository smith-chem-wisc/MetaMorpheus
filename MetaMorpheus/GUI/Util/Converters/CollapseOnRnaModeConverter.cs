using System;
using System.Globalization;
using System.Windows.Data;
using GuiFunctions;

namespace MetaMorpheusGUI;

public class CollapseOnRnaModeConverter : IValueConverter
{
    public object Convert(object value, Type targetType, object parameter, CultureInfo culture)
    {
        // Collapse the element if RNA mode is enabled
        return GuiGlobalParamsViewModel.Instance == null || !GuiGlobalParamsViewModel.Instance.IsRnaMode
            ? System.Windows.Visibility.Visible
            : System.Windows.Visibility.Collapsed;
    }
    public object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
    {
        throw new NotImplementedException();
    }
}