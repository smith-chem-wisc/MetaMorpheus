using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public enum PseudoMs2ConstructionType
    {
        MzPeak, //for analysis where MS2 peak tracing is performed at mz level
        Mass, //for analysis where MS2 peak tracing is performed at isotopic envelope level after deconvolution (each peak gets a mass and a charge)
    }
}
