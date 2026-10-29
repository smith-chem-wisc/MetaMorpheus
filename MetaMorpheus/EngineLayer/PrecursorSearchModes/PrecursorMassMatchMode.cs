namespace EngineLayer
{
    /// <summary>
    /// Selects which precursor mass is used to choose theoretical proteoform candidates during search.
    /// </summary>
    public enum PrecursorMassMatchMode
    {
        /// <summary>
        /// Match candidates by the monoisotopic precursor mass (the historical default). For high-mass
        /// proteoforms the monoisotopic peak is often undetectable, which can produce off-by-N errors.
        /// </summary>
        Monoisotopic,

        /// <summary>
        /// Match candidates by the most abundant observed isotopic peak (or, for isotopically
        /// unresolved high-mass species, the centroid/average mass). This is the experimentally most
        /// detectable peak in each charge envelope and avoids monoisotopic off-by-N errors.
        /// </summary>
        MostAbundant
    }
}
