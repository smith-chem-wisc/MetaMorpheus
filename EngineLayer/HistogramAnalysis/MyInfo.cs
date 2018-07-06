namespace EngineLayer.HistogramAnalysis
{
    internal class MyInfo
    {
        internal string Infostring;

        public MyInfo(double massShift, string infostring)
        {
            MassShift = massShift;
            Infostring = infostring;
        }

        public double MassShift { get; internal set; }
    }
}