using SpectralAveraging;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MetaMorpheusGUI
{
    public class SpectraFileProcessingTypeToNumberToAverageReadOnlyConverter : BaseValueConverter<SpectraFileProcessingTypeToNumberToAverageReadOnlyConverter>
    {
        public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            if (value.GetType() == typeof(SpectraFileAveragingType))
            {
                switch (value)
                {
                    case SpectraFileAveragingType.AverageAll:
                        return true;

                    case SpectraFileAveragingType.AverageEverynScans:
                        return false;

                    case SpectraFileAveragingType.AverageEverynScansWithOverlap:
                        return false;

                    case SpectraFileAveragingType.AverageDdaScans:
                        return false;

                    case SpectraFileAveragingType.AverageDdaScansWithOverlap:
                        return false;
                }
            }

            return false;
        }

        public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
        {
            throw new NotImplementedException();
        }
    }
}
