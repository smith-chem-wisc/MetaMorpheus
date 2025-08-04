using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class DIAparameters
    {
        public XicConstructor Ms1XicConstructor { get; set; }
        public XicConstructor Ms2XicConstructor { get; set; } 
        public PfGroupingEngine PfGroupingEngine { get; set; }
        public PseudoMs2ConstructionType PseudoMs2ConstructionType { get; set; }
        public bool CombineFragments { get; set; } 

        public DIAparameters(XicConstructor ms1XicConstructor, XicConstructor ms2XicConstructor, PfGroupingEngine pfGroupingEngine, PseudoMs2ConstructionType pseudoMs2ConstructionType, bool combineFragments = false)
        {
            Ms1XicConstructor = ms1XicConstructor;
            Ms2XicConstructor = ms2XicConstructor;
            PfGroupingEngine = pfGroupingEngine;
            PseudoMs2ConstructionType = pseudoMs2ConstructionType;
            CombineFragments = combineFragments;
        }
    }
}
