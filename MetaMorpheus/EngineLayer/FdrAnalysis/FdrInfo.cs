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
        /// True when <see cref="PEP_QValue"/> at this level comes from a PEP model trained for this match's digestion
        /// agent. Set by <see cref="FdrAnalysisEngine"/> for each digestion agent it trains or declines to train; false
        /// when PEP was not trained (too few matches, or PEP not requested) or training failed, in which case a PEP
        /// q-value may still have been computed, but only as a ranking of untrained PEPs.
        /// </summary>
        public bool PepTrained { get; set; }

        /// <summary>
        /// Creates a new FdrInfo object where Q-Values and PEP_Qvalues are set to 2 by default
        /// This is done to avoid situations where q-values aren't calcualted for a given peptides, but it is still
        /// reported in the final results. 
        /// </summary>
        public FdrInfo()
        {
            QValue = 2;
            PEP_QValue = 2;
            QValueNotch = 2;
        }

        /// <summary>
        /// A new FdrInfo holding the same values. Writing to the copy leaves this one unchanged. Every property is a
        /// number or a bool, so a member-wise copy is a full copy, and a property added later is copied without a change here.
        /// </summary>
        public FdrInfo Clone() => (FdrInfo)MemberwiseClone();
    }
}