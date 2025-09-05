using System;
using System.Globalization;

namespace MetaMorpheusGUI;

public class StringIsNullOrEmptyToBoolConverter : BaseValueConverter<StringIsNullOrEmptyToBoolConverter>
{
    public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
    {
        if (value is string str)
        {
            return string.IsNullOrEmpty(str);
        }
        return true; // or false, depending on how you want to handle non-string inputs
    }
    public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
    {
        throw new NotImplementedException();
    }
}