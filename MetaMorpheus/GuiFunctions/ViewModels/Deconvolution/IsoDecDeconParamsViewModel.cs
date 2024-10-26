using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;

namespace GuiFunctions
{
    public class IsoDecDeconParamsViewModel : DeconParamsViewModel
    {
        public override DeconvolutionParameters Parameters { get; protected set; }
    }
}
