namespace EngineLayer.SpectrumMatch
{
    /// <summary>
    /// The FDR values a match or peptide can be judged by. <see cref="CommonParameters.IdentificationFilter"/> is
    /// <see cref="PepQValue"/> or <see cref="QValueNotch"/>.
    /// </summary>
    public enum FilterType
    {
        /// <summary>
        /// The pooled q-value, counted over every notch together. Not an identification filter:
        /// <see cref="CommonParameters.IdentificationFilter"/> and <see cref="FilteredPsms"/> refuse it. It stays in the
        /// enum so that a task TOML naming it fails to load with a message that says why.
        /// </summary>
        QValue,

        /// <summary>
        /// The PEP q-value: matches ranked by their posterior error probability, decoys counted down that ranking. A match
        /// whose PEP was not trained is judged by its q-value notch instead.
        /// </summary>
        PepQValue,

        /// <summary>
        /// The q-value notch: the q-value counted within each mass-difference notch. The default identification filter.
        /// </summary>
        QValueNotch
    }
}
