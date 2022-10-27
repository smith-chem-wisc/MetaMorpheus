using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace MetaMorpheusGUI
{
    public class CheckboxUncheckConverter : BaseValueConverter<CheckboxUncheckConverter>
    {
        public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            if ((bool)value)
            {
                return false;
            }
            else
            {
                return (bool)value;
            }
        }

        public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
        {
            if (!(bool)value)
            {
                return false;
            }
            else
            {
                return (bool)parameter;
            }
        }
    }
}
