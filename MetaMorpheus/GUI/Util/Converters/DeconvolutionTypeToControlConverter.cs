using System;
using System.Globalization;
using GuiFunctions;
using MassSpectrometry;

namespace MetaMorpheusGUI
{
    public class DeconvolutionTypeToControlConverter : BaseValueConverter<DeconvolutionTypeToControlConverter>
    {
        public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            // Assuming value is of an enum type DeconvolutionType
            if (value is DeconParamsViewModel viewModel)
            {
                switch (viewModel.DeconvolutionType)
                {
                    case DeconvolutionType.ClassicDeconvolution:
                        return new ClassicDeconParamsControl() { DataContext = value as ClassicDeconParamsViewModel };
                
                    case DeconvolutionType.ExampleNewDeconvolutionTemplate:
                    default:
                        throw new ArgumentException("Invalid DeconvolutionType", nameof(value));
                }
            }
            throw new ArgumentException("Invalid value type", nameof(value));
        }

        public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
        {
            throw new NotImplementedException();
        }
    }
}
