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

        public double PEP { get; set; }
        public double PEP_QValue { get; set; }
    }
}