using MzLibUtil;
using System;
using System.Globalization;
using System.Windows;
using System.Windows.Data;
using TaskLayer;

namespace MetaMorpheusGUI;

/// <summary>
/// Converts an enum value to a boolean by comparing it to a parameter.
/// </summary>
public class EnumToBooleanConverter : IMultiValueConverter
{
    public object Convert(object[] values, Type targetType, object parameter, CultureInfo culture)
    {
        if (values.Length != 2 || values[0] == null || values[1] == null)
            return false;

        var type0 = values[0].GetType();
        var type1 = values[1].GetType();

        if (type0.IsEnum && type1.IsEnum && type0 == type1)
        {
            return values[0].Equals(values[1]);
        }

        return false;
    }

    public object[] ConvertBack(object value, Type[] targetTypes, object parameter, CultureInfo culture)
    {
        return new object[] { Binding.DoNothing, Binding.DoNothing };
    }
}

/// <summary>
/// Converts a MassDiffAcceptorType to a Visibility value.
/// Used to determine if the custom text window should be displayed based upon the selected mass diff acceptor type
/// </summary>
public class MassDifferenceAcceptorTypeToCustomTextBoxVisibilityConverter : BaseValueConverter<MassDifferenceAcceptorTypeToCustomTextBoxVisibilityConverter>
{
    public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
    {
        if (value is MassDiffAcceptorType type)
        {
            if (type == MassDiffAcceptorType.Custom)
            {
                return Visibility.Visible;
            }
        }
        return Visibility.Collapsed;
    }

    public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
    {
        throw new NotImplementedException();
    }
}