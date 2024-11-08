using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;

namespace GuiFunctions
{
    // Coming Soon
    [ExcludeFromCodeCoverage]
    public sealed class IsoDecDeconParamsViewModel : DeconParamsViewModel
    {
        public override DeconvolutionParameters Parameters { get; protected set; }


        public override string ToString() => "IsoDec";
    }
}
