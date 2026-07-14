namespace EngineLayer
{
    /// <summary>
    /// Chooses which observed precursor mass a search matches candidates against.
    /// <para>
    /// <see cref="Ms2ScanWithSpecificMass"/> and <see cref="SpectralMatch"/> store only OBSERVATIONS —
    /// the deconvoluted monoisotopic mass and the most-abundant (tallest) isotopologue mass. Neither
    /// property changes meaning with the search being run, so reading
    /// <see cref="Ms2ScanWithSpecificMass.PrecursorMass"/> always gives the monoisotopic mass and can
    /// never silently return something else. The knowledge of WHICH observation a given search keys off
    /// belongs to the <see cref="MassDiffAcceptor"/>, which is the object that defines the search type;
    /// these extensions are the one place that asks it.
    /// </para>
    /// </summary>
    public static class PrecursorMassExtensions
    {
        /// <summary>
        /// The observed precursor mass used to select theoretical candidates: the monoisotopic mass for a
        /// monoisotopic acceptor, the most-abundant isotopologue mass for a most-abundant acceptor. Falls
        /// back to the monoisotopic mass when no most-abundant peak was observed — no deconvoluted
        /// envelope (e.g. a scan-header precursor, or a neutral mass read from a pre-deconvoluted file).
        /// </summary>
        public static double GetPrecursorMassForSearch(this Ms2ScanWithSpecificMass scan, MassDiffAcceptor massDiffAcceptor)
            => Select(scan.PrecursorMass, scan.PrecursorMostAbundantMass, massDiffAcceptor);

        /// <inheritdoc cref="GetPrecursorMassForSearch(Ms2ScanWithSpecificMass, MassDiffAcceptor)"/>
        public static double GetPrecursorMassForSearch(this SpectralMatch psm, MassDiffAcceptor massDiffAcceptor)
            => Select(psm.ScanPrecursorMass, psm.ScanPrecursorMostAbundantMass, massDiffAcceptor);

        /// <summary>
        /// The observed monoisotopic mass implied by this match, with any isotope-level offset the
        /// acceptor allowed removed (see <see cref="MassDiffAcceptor.ToObservedMonoisotopicMass"/>).
        /// Subtracting the candidate's theoretical monoisotopic mass from this gives the CHEMICAL mass
        /// difference — what localization localizes, what G-PTM-D turns into a modification, and what
        /// calibration reads as instrument error. Equals <see cref="SpectralMatch.ScanPrecursorMass"/>
        /// in the default monoisotopic mode.
        /// </summary>
        public static double GetObservedMonoisotopicMass(this SpectralMatch psm, double theoreticalMonoisotopicMass, MassDiffAcceptor massDiffAcceptor)
        {
            // If the search matched on the monoisotopic observation — no acceptor, a monoisotopic acceptor,
            // or a most-abundant acceptor with no observed apex to fall back from — there is no isotope
            // offset to undo, and the monoisotopic mass is already the answer. Guarding here (rather than
            // handing the fallback mass to the acceptor) matters: ToObservedMonoisotopicMass assumes it is
            // given an apex mass and would otherwise subtract an apex offset from a monoisotopic one.
            if (massDiffAcceptor is not { MatchesMostAbundantPeak: true } || psm.ScanPrecursorMostAbundantMass <= 0)
            {
                return psm.ScanPrecursorMass;
            }

            return massDiffAcceptor.ToObservedMonoisotopicMass(psm.ScanPrecursorMostAbundantMass, theoreticalMonoisotopicMass);
        }

        private static double Select(double monoisotopicMass, double mostAbundantMass, MassDiffAcceptor massDiffAcceptor)
            => massDiffAcceptor is { MatchesMostAbundantPeak: true } && mostAbundantMass > 0
                ? mostAbundantMass
                : monoisotopicMass;
    }
}
