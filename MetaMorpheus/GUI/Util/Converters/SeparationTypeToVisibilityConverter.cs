using System;
using System.Globalization;
using System.Windows;

namespace MetaMorpheusGUI
{
    public class SeparationTypeToVisibilityConverter : BaseValueConverter<SeparationTypeToVisibilityConverter>
    {
        public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            return value is "CZE"
                ? Visibility.Collapsed
                : Visibility.Visible;
        }
        public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
        {
            throw new NotImplementedException();
        }
    }
}
