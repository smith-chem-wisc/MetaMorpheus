using System;
using System.Globalization;
using System.Windows.Data;

namespace MetaMorpheusGUI;

public class CollapseOnRnaModeConverter : IValueConverter
{
    public object Convert(object value, Type targetType, object parameter, CultureInfo culture)
    {
        // Collapse the element if RNA mode is enabled
        return UpdateGUISettings.Globals == null || !UpdateGUISettings.Globals.IsRnaMode
            ? System.Windows.Visibility.Visible
            : System.Windows.Visibility.Collapsed;
    }
    public object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
    {
        throw new NotImplementedException();
    }
}