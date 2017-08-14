using System.Collections.Generic;

namespace EngineLayer
{
    public static class ProductTypeToTerminusType
    {
        public static TerminusType IdentifyTerminusType(List<ProductType> lp)
        {
            if ((lp.Contains(ProductType.B) || lp.Contains(ProductType.BnoB1ions) || lp.Contains(ProductType.C) || lp.Contains(ProductType.Adot)) 
                && (lp.Contains(ProductType.Y) || lp.Contains(ProductType.Zdot) || lp.Contains(ProductType.X)))
            {
                return TerminusType.None;
            }
            else if (lp.Contains(ProductType.Y) || lp.Contains(ProductType.Zdot) || lp.Contains(ProductType.X))
            {
                return TerminusType.C;
            }
            else //if(lp.Contains(ProductType.B) || lp.Contains(ProductType.BnoB1ions) || lp.Contains(ProductType.C) || lp.Contains(ProductType.Adot))
            {
                return TerminusType.N;
            }
        }
    }
}
