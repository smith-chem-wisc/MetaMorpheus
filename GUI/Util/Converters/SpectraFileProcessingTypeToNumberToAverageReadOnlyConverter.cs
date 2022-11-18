using SpectralAveraging;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibSpectralAveraging;

namespace MetaMorpheusGUI
{
    public class SpectraFileProcessingTypeToNumberToAverageReadOnlyConverter : BaseValueConverter<SpectraFileProcessingTypeToNumberToAverageReadOnlyConverter>
    {
        public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            if (value.GetType() == typeof(SpectraFileProcessingType))
            {
                switch (value)
                {
                    case SpectraFileProcessingType.AverageAll:
                        return true;

                    case SpectraFileProcessingType.AverageEverynScans:
                        return false;

                    case SpectraFileProcessingType.AverageEverynScansWithOverlap:
                        return false;

                    case SpectraFileProcessingType.AverageDDAScans:
                        return false;

                    case SpectraFileProcessingType.AverageDDAScansWithOverlap:
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
