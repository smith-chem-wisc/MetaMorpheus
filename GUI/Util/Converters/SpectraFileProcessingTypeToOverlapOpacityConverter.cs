using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SpectralAveraging;

namespace MetaMorpheusGUI
{
    public class SpectraFileProcessingTypeToOverlapOpacityConverter : BaseValueConverter<SpectraFileProcessingTypeToOverlapOpacityConverter>
    {
        private double doNotDisplay = 0.5;

        public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            if (value.GetType() == typeof(SpectraFileAveragingType))
            {
                switch (value)
                {
                    case SpectraFileAveragingType.AverageAll:
                        return doNotDisplay;

                    case SpectraFileAveragingType.AverageEverynScans:
                        return doNotDisplay;

                    case SpectraFileAveragingType.AverageEverynScansWithOverlap:
                        return 1;

                    case SpectraFileAveragingType.AverageDdaScans:
                        return doNotDisplay;

                    case SpectraFileAveragingType.AverageDdaScansWithOverlap:
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
