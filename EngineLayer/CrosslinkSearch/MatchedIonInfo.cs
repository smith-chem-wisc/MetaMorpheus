namespace EngineLayer.CrosslinkSearch
{
    public class MatchedIonInfo
    {
        #region Public Constructors

        public MatchedIonInfo(int length)
        {
            MatchedIonMz = new double[length];
            MatchedIonIntensity = new double[length];
            MatchedIonName = new string[length];
        }

        #endregion Public Constructors

        #region Public Properties

        public double[] MatchedIonMz { get; set; }
        public string[] MatchedIonName { get; set; }
        public double[] MatchedIonIntensity { get; set; }

        #endregion Public Properties
    }
}