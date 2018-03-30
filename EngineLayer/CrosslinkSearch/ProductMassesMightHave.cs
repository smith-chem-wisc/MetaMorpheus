namespace EngineLayer.CrosslinkSearch
{
    public class ProductMassesMightHave
    {
        #region Public Constructors

        public ProductMassesMightHave(int length)
        {
            ProductMz = new double[length];
            ProductName = new string[length];
        }

        public ProductMassesMightHave()
        {
        }

        #endregion Public Constructors

        #region Public Properties

        public double[] ProductMz { get; set; }
        public string[] ProductName { get; set; }
        public int XlPos { get; set; }
        public int XlPos2 { get; set; }

        #endregion Public Properties
    }
}