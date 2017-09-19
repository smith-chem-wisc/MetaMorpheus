using System.Collections.Generic;

namespace EngineLayer
{
    public enum ProductType
    {
        Adot,
        B,
        BnoB1ions,
        C,
        X,
        Y,
        Zdot,
        None
    }

    public static class ProductTypeMethod
    {
        #region Public Methods

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

        public static List<List<ProductType>> SeparateIonsByTerminus(List<ProductType> ionTypes)
        {
            List<ProductType> nIons = new List<ProductType>();
            List<ProductType> cIons = new List<ProductType>();
            foreach (ProductType productType in ionTypes)
            {
                if (productType == ProductType.B || productType == ProductType.C)
                    nIons.Add(productType);
                else
                    cIons.Add(productType);
            }

            return new List<List<ProductType>> { nIons, cIons };
        }

        #endregion Public Methods
    }
}