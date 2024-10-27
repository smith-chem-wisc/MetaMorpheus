using System;
using System.Globalization;
using GuiFunctions;
using MassSpectrometry;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Converts a DeconvolutionType to a corresponding control.
    /// </summary>
    public class DeconvolutionTypeToControlConverter : BaseValueConverter<DeconvolutionTypeToControlConverter>
    {
        /// <summary>
        /// Converts a DeconvolutionType to a corresponding control.
        /// </summary>
        /// <param name="value">The value to convert, expected to be of type DeconParamsViewModel.</param>
        /// <param name="targetType">The type of the target property.</param>
        /// <param name="parameter">Optional parameter to be used in the converter logic.</param>
        /// <param name="culture">The culture to use in the converter.</param>
        /// <returns>A control corresponding to the DeconvolutionType.</returns>
        /// <exception cref="ArgumentException">Thrown when the value is not of type DeconParamsViewModel or the DeconvolutionType is invalid.</exception>
        public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            // Assuming value is of an enum type DeconvolutionType
            if (value is not DeconParamsViewModel viewModel)
                throw new ArgumentException("Invalid value type", nameof(value));

            switch (viewModel.DeconvolutionType)
            {
                case DeconvolutionType.ClassicDeconvolution:
                    return new ClassicDeconParamsControl() { DataContext = value as ClassicDeconParamsViewModel };

                case DeconvolutionType.ExampleNewDeconvolutionTemplate:
                default:
                    throw new ArgumentException("Invalid DeconvolutionType", nameof(value));
            }
        }

        /// <summary>
        /// Converts a value back to its source type. Not implemented.
        /// </summary>
        /// <param name="value">The value produced by the binding target.</param>
        /// <param name="targetType">The type to convert to.</param>
        /// <param name="parameter">Optional parameter to be used in the converter logic.</param>
        /// <param name="culture">The culture to use in the converter.</param>
        /// <returns>Throws NotImplementedException.</returns>
        /// <exception cref="NotImplementedException">Always thrown as this method is not implemented.</exception>
        public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
        {
            throw new NotImplementedException();
        }
    }
}
