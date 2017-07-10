namespace EngineLayer.CrosslinkSearch
{
    public class MatchedIonInfo
    {
        public double[] MatchedIonMz { get; set; }
        public string[] MatchedIonName { get; set; }
        public double[] MatchedIonIntensity { get; set; }

        public MatchedIonInfo(int length)
        {
            MatchedIonMz = new double[length];
            MatchedIonIntensity = new double[length];
            MatchedIonName = new string[length];
        }
    }
}