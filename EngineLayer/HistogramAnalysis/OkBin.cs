namespace EngineLayer.HistogramAnalysis
{
    internal class OkBin
    {
        internal double massShift;
        internal double sigma;
        internal int p;

        public OkBin(double massShift, double sigma, int p)
        {
            this.massShift = massShift;
            this.sigma = sigma;
            this.p = p;
        }
    }
}