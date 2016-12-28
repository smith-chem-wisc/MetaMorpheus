namespace IndexSearchAndAnalyze
{
    internal class OkBin
    {
        internal double massShift;
        internal double sigma;
        internal int p;

        public OkBin(double v1, double v2, int v3)
        {
            this.massShift = v1;
            this.sigma = v2;
            this.p = v3;
        }
    }
}