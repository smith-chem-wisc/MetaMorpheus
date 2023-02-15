using SpectralAveraging;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MetaMorpheusGUI
{
    public class OutlierRejectionTypeToPercentileOpacityConverter : BaseValueConverter<OutlierRejectionTypeToPercentileOpacityConverter>
    {
        private double doNotDisplay = 0.5;
        public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            if (value.GetType() == typeof(OutlierRejectionType))
            {
                if ((OutlierRejectionType)value == OutlierRejectionType.PercentileClipping)
                    return 1;
                else
                    return doNotDisplay;
            }

            return 1;
        }

        public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
        {
            throw new NotImplementedException();
        }
    }
}
