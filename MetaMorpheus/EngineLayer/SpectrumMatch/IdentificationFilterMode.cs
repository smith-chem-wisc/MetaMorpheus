namespace EngineLayer.SpectrumMatch
{
    /// <summary>
    /// Which FDR value decides whether a spectral match or peptide counts as identified, everywhere a search uses
    /// that decision: the PSM and peptide tables, protein inference and protein FDR, quantification, the pruned
    /// database, the spectral library, the mzIdentML and pepXML files, and the summaries.
    /// </summary>
    /// <remarks>
    /// A named filter (anything but <see cref="Thresholds"/>) compares strictly (value &lt; threshold). PEP is trained
    /// per digestion agent, and in a semi- or non-specific search per cleavage-specificity category. Where it was not trained, both PEP filters fall back
    /// for those matches only, and the search says so in results.txt, naming the digestion agent and the category,
    /// and warns.
    /// </remarks>
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
        /// The PEP q-value, strictly below the PEP q-value threshold. Where PEP was not trained, the q-value instead,
        /// read as MetaMorpheus's q-value filter reads it (the q-value and the q-value notch), below the same
        /// threshold.
        /// </summary>
        PepQValue,

        /// <summary>
        /// The PEP q-value where PEP was trained, otherwise the q-value notch alone, strictly below the smaller of the
        /// two thresholds. It differs from <see cref="PepQValue"/> in two ways. An untrained digestion agent falls back
        /// to the notch alone, which a match in a sparsely populated mass-difference notch can pass while its q-value
        /// fails. And the threshold is the smaller of the two, so it works whichever threshold box the GUI fills in.
        /// </summary>
        PepQValueElseQValueNotch
    }
}
