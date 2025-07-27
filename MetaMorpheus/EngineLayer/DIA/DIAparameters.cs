using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class DIAparameters
    {
        public DIAAnalysisType DIAAnalysisType { get; set; } 
        public XicConstructor XicConstructor { get; set; }
        public PfGroupingEngine PfGroupingEngine { get; set; }
        public PseudoMs2ConstructionType PseudoMs2ConstructionType { get; set; }

        public DIAparameters(DIAAnalysisType diaAnalysisType, XicConstructor xicConstructor, PfGroupingEngine pfGroupingEngine, PseudoMs2ConstructionType pseudoMs2ConstructionType)
        {
            DIAAnalysisType = diaAnalysisType;
            XicConstructor = xicConstructor;
            PfGroupingEngine = pfGroupingEngine;
            PseudoMs2ConstructionType = pseudoMs2ConstructionType;
        }
    }
}
