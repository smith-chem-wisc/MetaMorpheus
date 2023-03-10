namespace EngineLayer.HistogramAnalysis
{
    internal class OkBin
    {
        internal double MassShift;
        internal double Sigma;
        internal int P;

        public OkBin(double massShift, double sigma, int p)
        {
            MassShift = massShift;
            Sigma = sigma;
            P = p;
        }
    }
}