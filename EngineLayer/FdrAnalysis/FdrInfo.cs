namespace EngineLayer.FdrAnalysis
{
    public class FdrInfo
    {
        public int CumulativeTarget { get; set; }
        public int CumulativeDecoy { get; set; }
        public int CumulativeTargetNotch { get; set; }
        public int CumulativeDecoyNotch { get; set; }
        public double QValue { get; set; }
        public double QValueNotch { get; set; }
        public bool CalculateEValue { get; set; }
        public double MaximumLikelihood { get; set; }
        public double EValue { get; set; }
        public double EScore { get; set; }
    }
}