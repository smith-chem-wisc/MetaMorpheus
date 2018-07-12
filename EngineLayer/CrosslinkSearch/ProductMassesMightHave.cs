namespace EngineLayer.CrosslinkSearch
{
    public class ProductMassesMightHave
    {
        public ProductMassesMightHave(int length)
        {
            ProductMz = new double[length];
            ProductName = new string[length];
        }

        public ProductMassesMightHave()
        {
        }

        public double[] ProductMz { get; set; }
        public string[] ProductName { get; set; }
        public int XlPos { get; set; }
        public int XlPos2 { get; set; }
    }
}