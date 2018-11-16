using Proteomics.Fragmentation;
using System.Collections.Generic;
namespace MetaMorpheusGUI
{
    class CustomFragmentation
    {
        public bool BIons { get; set; }
        public bool CIons { get; set; }
        public bool YIons { get; set; }
        public bool ZDotIons { get; set; }

        public static List<ProductType> CustomFragmentationIons(CustomFragmentation c)
        {
            var fragmentationIons = new List<ProductType>();

            if (c.BIons)
            {
                fragmentationIons.Add(ProductType.b);
            }
            if (c.CIons)
            {
                fragmentationIons.Add(ProductType.c);
            }
            if (c.YIons)
            {
                fragmentationIons.Add(ProductType.y);
            }
            if (c.ZDotIons)
            {
                fragmentationIons.Add(ProductType.zDot);
            }

            return fragmentationIons;
        }
    }
}
