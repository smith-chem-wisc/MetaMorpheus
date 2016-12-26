namespace FragmentGeneration
{
    internal class MyInfo
    {
        internal string infostring;

        public MyInfo(double v1, string v2)
        {
            this.MassShift = v1;
            this.infostring = v2;
        }

        public double MassShift { get; internal set; }
    }
}