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
        PepQValue
    }
}
