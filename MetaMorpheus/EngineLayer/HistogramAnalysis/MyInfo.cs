namespace EngineLayer.HistogramAnalysis
{
    internal class MyInfo
    {
        internal string infostring;

        public MyInfo(double MassShift, string infostring)
        {
            this.MassShift = MassShift;
            this.infostring = infostring;
        }

        public double MassShift { get; internal set; }
    }
}