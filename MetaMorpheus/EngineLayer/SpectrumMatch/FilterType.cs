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
        /// Used only when a caller asks for it by name; MetaMorpheus never infers it.
        /// </summary>
        QValueNotch
    }
}
