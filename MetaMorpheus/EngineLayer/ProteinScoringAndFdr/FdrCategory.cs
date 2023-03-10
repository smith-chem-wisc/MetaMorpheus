using System;
using System.Collections.Generic;
using System.Text;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer
{
    public enum FdrCategory
    {
        //Cleavage Specificity
        FullySpecific = 0,
        SemiSpecific = 1,
        NonSpecific = 2,

        //New category here
    }
}
