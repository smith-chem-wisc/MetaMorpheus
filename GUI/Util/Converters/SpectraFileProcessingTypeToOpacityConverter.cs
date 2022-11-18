using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibSpectralAveraging;
using SpectralAveraging;

namespace MetaMorpheusGUI
{
    public class SpectraFileProcessingTypeToOpacityConverter : BaseValueConverter<SpectraFileProcessingTypeToOpacityConverter>
    {
        private double doNotDisplay = 0.5;

        public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            if (value.GetType() == typeof(SpectraFileProcessingType))
            {
                switch (value)
                {
                    case SpectraFileProcessingType.AverageAll:
                        return doNotDisplay;

                    case SpectraFileProcessingType.AverageEverynScans:
                        return 1;

                    case SpectraFileProcessingType.AverageEverynScansWithOverlap:
                        return 1;

                    case SpectraFileProcessingType.AverageDDAScans:
                        return 1;

                    case SpectraFileProcessingType.AverageDDAScansWithOverlap:
                        return 1;
                }
            }

            return 1;
        }

        public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
        {
            throw new NotImplementedException();
        }
    }
}
