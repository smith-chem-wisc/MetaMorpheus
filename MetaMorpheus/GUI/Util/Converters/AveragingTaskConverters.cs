using SpectralAveraging;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MetaMorpheusGUI
{
    

    public class SpectraFileAveragingTypeToOpacityConverter : BaseValueConverter<SpectraFileAveragingTypeToOpacityConverter>
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
                        return 1;

                    case SpectraFileAveragingType.AverageEverynScansWithOverlap:
                        return 1;

                    case SpectraFileAveragingType.AverageDdaScans:
                        return 1;

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

    public class SpectraFileAveragingTypeToNumberToAverageReadOnlyConverter : BaseValueConverter<SpectraFileAveragingTypeToNumberToAverageReadOnlyConverter>
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

    public class SpectraFileAveragingTypeToNumberOpacityConverter : BaseValueConverter<SpectraFileAveragingTypeToNumberOpacityConverter>
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
                        return 1;

                    case SpectraFileAveragingType.AverageEverynScansWithOverlap:
                        return 1;

                    case SpectraFileAveragingType.AverageDdaScans:
                        return 1;

                    case SpectraFileAveragingType.AverageDdaScansWithOverlap:
                        return 1;
                }
            }

            return false;
        }

        public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
        {
            throw new NotImplementedException();
        }
    }

    public class OutlierRejectionTypeToSigmaReadOnlyConverter : BaseValueConverter<OutlierRejectionTypeToSigmaReadOnlyConverter>
    {
        public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            if (value.GetType() == typeof(OutlierRejectionType))
            {
                OutlierRejectionType val = (OutlierRejectionType)value;
                if (val is OutlierRejectionType.SigmaClipping or OutlierRejectionType.AveragedSigmaClipping or OutlierRejectionType.WinsorizedSigmaClipping)
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

    public class OutlierRejectionTypeToPercentileReadOnlyConverter : BaseValueConverter<OutlierRejectionTypeToPercentileReadOnlyConverter>
    {
        public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            if (value.GetType() == typeof(OutlierRejectionType))
            {
                if ((OutlierRejectionType)value == OutlierRejectionType.PercentileClipping)
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

    public class DataStructureToStringConverter : BaseValueConverter<DataStructureToStringConverter>
    {
        public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            string val = value.ToString();

            // spectra file processing specific
            if (value.GetType() == typeof(SpectraFileAveragingType))
            {
                if (val == "AverageEverynScansWithOverlap")
                    return "Average Every n Scans";
                else if (val == "AverageDdaScansWithOverlap")
                    return "Average DDA Scans";
            }

            // output type specific
            if (value.GetType() == typeof(OutputType))
            {
                if (val == "txt")
                    return "Text File";
                else if (val == "mzML")
                    return "MzML file";
            }


            // empty
            if (string.IsNullOrWhiteSpace(val))
                return "";

            // split all others at the capital letter
            StringBuilder newText = new StringBuilder(val.Length * 2);
            newText.Append(val[0]);
            for (int i = 1; i < val.Length; i++)
            {
                if (char.IsUpper(val[i]) && val[i - 1] != ' ')
                    newText.Append(' ');
                newText.Append(val[i]);
            }
            return new String(newText.ToString());

        }

        public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
        {
            throw new NotImplementedException();
        }
    }
}
