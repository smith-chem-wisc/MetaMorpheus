namespace EngineLayer.SpectrumMatch
{
    /// <summary>
    /// The FDR values a match or peptide can be judged by: the three choices for
    /// <see cref="CommonParameters.IdentificationFilter"/>, each held to <see cref="CommonParameters.FilterThreshold"/>.
    /// </summary>
    public enum FilterType
    {
        /// <summary>
        /// The pooled q-value alone, counted over every mass-difference notch together. Before the identification filter
        /// was named, a "q-value" filter meant the pooled q-value and the q-value notch together; this is the pooled
        /// q-value alone.
        /// </summary>
        QValue,

        /// <summary>
        /// The PEP q-value: matches ranked by their posterior error probability, decoys counted down that ranking. A match
        /// whose PEP was not trained is judged by its q-value notch instead. The default bottom-up
        /// (<see cref="CommonParameters.DefaultIdentificationFilter"/>).
        /// </summary>
        PepQValue,

        /// <summary>
        /// The q-value notch: the q-value counted within each mass-difference notch. The default for top-down and RNA
        /// (<see cref="CommonParameters.DefaultIdentificationFilter"/>).
        /// </summary>
        QValueNotch,

        /// <summary>A pooled q-value and q-value-notch filter, both held to the same threshold.</summary>
        QValueAndQValueNotch
    }
}
