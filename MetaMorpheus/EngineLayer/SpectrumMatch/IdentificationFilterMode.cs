namespace EngineLayer.SpectrumMatch
{
    /// <summary>
    /// Which FDR value decides whether a spectral match or peptide counts as identified, everywhere a search uses
    /// that decision: the PSM and peptide tables, protein inference and protein FDR, quantification, the pruned
    /// database, the spectral library, the mzIdentML and pepXML files, and the summaries.
    /// </summary>
    public enum IdentificationFilterMode
    {
        /// <summary>
        /// The filter is inferred from the thresholds, as MetaMorpheus always has: the PEP q-value when the PEP
        /// q-value threshold is below the q-value threshold, otherwise the q-value and the q-value notch, both
        /// inclusive (value &lt;= threshold). The default.
        /// </summary>
        Thresholds,

        /// <summary>
        /// The q-value notch alone, strictly below the q-value threshold.
        /// </summary>
        QValueNotch,

        /// <summary>
        /// The PEP q-value, strictly below the PEP q-value threshold. A search in which PEP was not trained for a
        /// digestion agent fails, naming it.
        /// </summary>
        PepQValue,

        /// <summary>
        /// The PEP q-value where PEP was trained for the match's digestion agent, otherwise the q-value notch; strictly
        /// below the smaller of the two thresholds. Every fallback is reported.
        /// </summary>
        PepQValueElseQValueNotch
    }
}
