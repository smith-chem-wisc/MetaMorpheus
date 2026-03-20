using System;
using System.Globalization;

namespace MetaMorpheusGUI;

public class IntegerGreaterThanZeroToBoolConverter : BaseValueConverter<IntegerGreaterThanZeroToBoolConverter>
{
    public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
    {
        if (value is int intValue)
        {
            return intValue > 0;
        }
        return false;
    }

    public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
    {
        throw new NotImplementedException();
    }
}
