namespace EngineLayer.GlycoSearch
{
    /// <summary>
    /// Frozen positions inside <see cref="Glycan.AllOxoniumIons"/> that legacy production code
    /// indexes by number. Naming these is the only refactor; the ordering of AllOxoniumIons
    /// itself MUST NOT change. The <c>Reserved_Indices_Match_AllOxoniumIons_Positions</c> test
    /// guards the ordering by asserting the scaled-int mass at each index below.
    ///
    /// Ion assignments are well-established in the glycoproteomics literature
    /// (see Halim et al., Anal Chem 2014 for oxonium fragmentation profiles;
    /// Polasky et al., J Proteome Res 2024 for experimentally validated diagnostic ion sets;
    /// Toghi Eshghi et al., Anal Chem 2016 for the 138/144 ratio as an N- vs O-glycopeptide classifier).
    /// </summary>
    internal static class OxoniumIonReservedIndices
    {
        /// <summary>HexNAc fragment at 138.055 (numerator of the 138/144 ratio used to
        /// discriminate N-glycan from O-glycan; both 138 and 144 are HexNAc-derived, the
        /// <i>ratio</i> is diagnostic).</summary>
        public const int R138 = 4;

        /// <summary>HexNAc fragment at 144.066 (denominator of the same ratio).</summary>
        public const int R144 = 5;

        /// <summary>HexNAc oxonium 204.087 — canonical glycopeptide marker; spectra lacking
        /// this peak are filtered out in FindNGlycan and MatchGlycopeptide.</summary>
        public const int HexNAc204 = 9;

        /// <summary>NeuAc-H2O 274.093 (sialylation marker).</summary>
        public const int NeuAc274 = 10;

        /// <summary>NeuAc 292.103 (sialic acid marker).</summary>
        public const int NeuAc292 = 12;

        /// <summary>HexHexNAc 366.140 (LacNAc-type composite diagnostic).</summary>
        public const int HexHexNAc366 = 14;
    }
}
