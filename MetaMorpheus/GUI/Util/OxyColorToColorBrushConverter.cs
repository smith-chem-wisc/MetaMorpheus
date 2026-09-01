using System;
using System.Globalization;
using System.Windows.Media;
using OxyPlot;
using OxyPlot.Wpf;

namespace MetaMorpheusGUI;

public class OxyColorToColorBrushConverter : BaseValueConverter<OxyColorToColorBrushConverter>
{
    public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
    {
        if (object.ReferenceEquals(value, null))
            return null;
        if (value is OxyColor oxyColor)
            return oxyColor.ToBrush();
        return null;
    }

    public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
    {
        if (object.ReferenceEquals(value, null))
            return null;
        if (value is Brush brush)
            return brush.ToOxyColor();
        return null;
    }
}