using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SpectralAveraging;

namespace MetaMorpheusGUI
{
    public class OutlierRejectionTypeToPercentileReadOnlyConverter : BaseValueConverter<OutlierRejectionTypeToPercentileReadOnlyConverter>
    {
        public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            if (value.GetType() == typeof(RejectionType))
            {
                if ((RejectionType)value == RejectionType.PercentileClipping)
                    return false;
                else 
                    return true;
            }

            return 1;
        }

        public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
        {
            throw new NotImplementedException();
        }
    }
}
