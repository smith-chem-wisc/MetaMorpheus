using System;
using System.Globalization;
using System.IO;

namespace MetaMorpheus
{
    public class ProductCaps
    {
        private static readonly double[,] PRODUCT_CAP_MASSES = new double[Enum.GetNames(typeof(ProductType)).Length, Enum.GetNames(typeof(MassType)).Length];

        private static readonly ProductCaps instance = new ProductCaps();

        private ProductCaps()
        {
            using (StreamReader product_caps = new StreamReader(Path.Combine(Path.GetDirectoryName(Environment.GetCommandLineArgs()[0]), "product_caps.tsv")))
            {
                product_caps.ReadLine();

                while (product_caps.Peek() != -1)
                {
                    string line = product_caps.ReadLine();
                    string[] fields = line.Split('\t');

                    ProductType product_type = (ProductType)Enum.Parse(typeof(ProductType), fields[0], true);
                    double cap_monoisotopic_mass = double.Parse(fields[1], CultureInfo.InvariantCulture);
                    PRODUCT_CAP_MASSES[(int)product_type, (int)MassType.Monoisotopic] = cap_monoisotopic_mass;
                }
            }
        }

        public static ProductCaps Instance
        {
            get { return instance; }
        }

        public double this[ProductType productType, MassType massType]
        {
            get { return PRODUCT_CAP_MASSES[(int)productType, (int)massType]; }
        }
    }
}