using System.Collections.Generic;

namespace EngineLayer
{
    static class ProductTypeToTerminusType
    {
        public static TerminusType IdentifyTerminusType(List<ProductType> lp)
        {
            if ((lp.Contains(ProductType.B) || lp.Contains(ProductType.BnoB1ions) || lp.Contains(ProductType.C)) && (lp.Contains(ProductType.Y) || lp.Contains(ProductType.Zdot)))
            {
                return TerminusType.None;
            }
            else if (lp.Contains(ProductType.Y) || lp.Contains(ProductType.Zdot))
            {
                return TerminusType.C;
            }
            else //if(lp.Contains(ProductType.B) || lp.Contains(ProductType.BnoB1ions) || lp.Contains(ProductType.C))
            {
                return TerminusType.N;
            }
        }
    }
}
