using System;
using MassSpectrometry;

namespace GuiFunctions
{
    /// <summary>
    /// Provides additional functionality for the MzLib library.
    /// </summary>
    public static class MzLibExtensions
    {

        /// <summary>
        /// Converts the given <see cref="DeconvolutionParameters"/> to a <see cref="DeconParamsViewModel"/>.
        /// </summary>
        /// <param name="parameters">The deconvolution parameters to convert.</param>
        /// <returns>A <see cref="DeconParamsViewModel"/> representing the given parameters.</returns>
        /// <exception cref="NotImplementedException">
        /// Thrown when the type of <paramref name="parameters"/> is not supported.
        /// </exception>
        public static DeconParamsViewModel ToViewModel(this DeconvolutionParameters parameters)
        {
            if (parameters is ClassicDeconvolutionParameters classicParams)
            {
                return new ClassicDeconParamsViewModel(classicParams);
            }
            else if (parameters is IsoDecDeconvolutionParameters isoParams)
            {
                return new IsoDecDeconParamsViewModel(isoParams);
            }
            else
            {
                throw new NotImplementedException();
            }
        }
    }
}
