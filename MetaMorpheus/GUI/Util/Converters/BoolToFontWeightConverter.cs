using System;
using System.Globalization;
using System.Windows;

namespace MetaMorpheusGUI;

public class BoolToFontWeightConverter : BaseValueConverter<BoolToFontWeightConverter>
{
    public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
    {
        if (value is bool boolValue && boolValue)
        {
            return FontWeights.Bold;
        }
        return FontWeights.Normal;
    }

    public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
    {
        throw new NotImplementedException();
    }
}
