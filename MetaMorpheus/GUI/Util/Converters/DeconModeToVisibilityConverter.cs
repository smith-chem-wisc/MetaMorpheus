using System;
using System.Globalization;
using System.Windows;
using GuiFunctions;

namespace MetaMorpheusGUI;

public class DeconModeToVisibilityConverter : BaseValueConverter<DeconModeToVisibilityConverter>
{
    public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
    {
        return value is DeconvolutionMode.FullSpectrum 
            ? Visibility.Visible 
            : Visibility.Collapsed;
    }

    public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
    {
        throw new NotImplementedException();
    }
}