using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;

namespace GuiFunctions
{
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
            //else if (parameters is IsoDeconvolutionParameters isoParams)
            //{
            //    return new IsoDecDeconParamsViewModel(isoParams);
            //}
            else
            {
                throw new NotImplementedException();
            }
        }
    }
}
