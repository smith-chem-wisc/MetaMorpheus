using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer
{
    public static class FdrClassifier
    {
        public static FdrCategory GetCleavageSpecificityCategory(CleavageSpecificity cleavageSpecificity)
        {
            if (cleavageSpecificity == CleavageSpecificity.Full)
            {
                return FdrCategory.FullySpecific;
            }
            else if (cleavageSpecificity == CleavageSpecificity.Semi)
            {
                return FdrCategory.SemiSpecific;
            }
            else if (cleavageSpecificity == CleavageSpecificity.None)
            {
                return FdrCategory.NonSpecific;
            }
            else
            {
                throw new NotImplementedException("Cleavage specificity '" + cleavageSpecificity + "' has not been immplemented for local FDR calculations.");
            }
        }
    }
}
