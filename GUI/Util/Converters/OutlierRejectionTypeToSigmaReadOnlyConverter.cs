using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SpectralAveraging;

namespace MetaMorpheusGUI
{
    public class OutlierRejectionTypeToSigmaReadOnlyConverter : BaseValueConverter<OutlierRejectionTypeToSigmaReadOnlyConverter>
    {
        public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            if (value.GetType() == typeof(RejectionType))
            {
                RejectionType val = (RejectionType)value;
                if (val is RejectionType.SigmaClipping or RejectionType.AveragedSigmaClipping or RejectionType.WinsorizedSigmaClipping)
                    return false;
                else
                    return true;
            }

            return false;
        }

        public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
        {
            throw new NotImplementedException();
        }
    }
}
