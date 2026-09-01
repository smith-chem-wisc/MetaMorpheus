using System;
using System.Globalization;
using System.Windows.Data;

namespace MetaMorpheusGUI;

public class SpectrumDescriptorToolTipConverter : IValueConverter
{
    public object Convert(object value, Type targetType, object parameter, CultureInfo culture)
    {
        var valueString = value?.ToString();
        if (valueString == null)
            return null;

        if (valueString.Contains("Ambiguity"))
            return "J. Proteome Res. 2018, 17, 3, 1321–1325";

        return null;
    }
    public object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture) => throw new NotImplementedException();
}