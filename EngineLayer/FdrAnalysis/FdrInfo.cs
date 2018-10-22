namespace EngineLayer.FdrAnalysis
{
    public class FdrInfo
    {
        public double CumulativeTarget { get; set; }
        public double CumulativeDecoy { get; set; }
        public double CumulativeTargetNotch { get; set; }
        public double CumulativeDecoyNotch { get; set; }
        public double QValue { get; set; }
        public double QValueNotch { get; set; }
        public bool CalculateEValue { get; set; }
        public double MaximumLikelihood { get; set; }
        public double EValue { get; set; }
        public double EScore { get; set; }
    }
}