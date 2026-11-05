namespace EngineLayer.SpectrumMatch
{
    /// <summary>
    /// Specifies the type of filtering used for PSM/peptide/protein FDR calculations.
    /// </summary>
    public enum FilterType
    {
        /// <summary>
        /// Filter based on q-value (classic FDR)
        /// </summary>
        QValue,

        /// <summary>
        /// Filter based on PEP q-value (posterior error probability)
        /// </summary>
        PepQValue,

        /// <summary>
        /// Filter based on the q-value notch alone (the q-value computed within each mass-difference notch).
        /// Used only when a caller or the search's identification filter (<see cref="IdentificationFilterMode"/>) names it;
        /// MetaMorpheus never infers it from the thresholds.
        /// </summary>
        QValueNotch
    }
}
