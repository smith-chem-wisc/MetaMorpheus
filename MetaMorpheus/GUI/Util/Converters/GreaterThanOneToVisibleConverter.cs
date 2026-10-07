using System;
using System.Globalization;

namespace MetaMorpheusGUI;

public class GreaterThanOneToVisibleConverter : BaseValueConverter<GreaterThanOneToVisibleConverter>
{
    public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
    {
        if (value is int intValue)
        {
            return intValue > 1 ? System.Windows.Visibility.Visible : System.Windows.Visibility.Collapsed;
        }
        return System.Windows.Visibility.Collapsed;
    }

    public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
    {
        // Not typically used; return Binding.DoNothing to indicate no conversion back
        return System.Windows.Data.Binding.DoNothing;
    }
}