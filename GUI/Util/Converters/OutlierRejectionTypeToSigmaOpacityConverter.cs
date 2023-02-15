using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SpectralAveraging;

namespace MetaMorpheusGUI
{
    public class OutlierRejectionTypeToSigmaOpacityConverter : BaseValueConverter<OutlierRejectionTypeToSigmaOpacityConverter>
    {
        private double doNotDisplay = 0.5;
        public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            if (value.GetType() == typeof(OutlierRejectionType))
            {
                OutlierRejectionType val = (OutlierRejectionType)value;
                if (val is OutlierRejectionType.SigmaClipping or OutlierRejectionType.AveragedSigmaClipping or OutlierRejectionType.WinsorizedSigmaClipping)
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
