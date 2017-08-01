namespace EngineLayer.CrosslinkSearch
{
    public class ProductMassesMightHave
    {
        public double[] ProductMz { get; set; }
        public string[] ProductName { get; set; }

        public ProductMassesMightHave(int length)
        {
            ProductMz = new double[length];
            ProductName = new string[length];
        }

        public ProductMassesMightHave()
        {
        }
    }
}
