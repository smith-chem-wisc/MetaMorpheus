namespace EngineLayer.HistogramAnalysis
{
    internal class OkBin
    {
        #region Internal Fields

        internal double massShift;
        internal double sigma;
        internal int p;

        #endregion Internal Fields

        #region Public Constructors

        public OkBin(double massShift, double sigma, int p)
        {
            this.massShift = massShift;
            this.sigma = sigma;
            this.p = p;
        }

        #endregion Public Constructors
    }
}