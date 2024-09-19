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

        /// <summary>
        /// Creates a new FdrInfo object where Q-Values and PEP_Qvalues are set to 2 by default
        /// This is done to avoid situations where q-values aren't calcualted for a given peptides, but it is still
        /// reported in the final results. 
        /// </summary>
        public FdrInfo()
        {
            QValue = 2;
            PEP_QValue = 2;
        }
    }
}