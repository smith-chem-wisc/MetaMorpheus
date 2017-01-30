namespace EngineLayer.Analysis
{
    internal class MyInfo
    {

        #region Internal Fields

        internal string infostring;

        #endregion Internal Fields

        #region Public Constructors

        public MyInfo(double MassShift, string infostring)
        {
            this.MassShift = MassShift;
            this.infostring = infostring;
        }

        #endregion Public Constructors

        #region Public Properties

        public double MassShift { get; internal set; }

        #endregion Public Properties

    }
}