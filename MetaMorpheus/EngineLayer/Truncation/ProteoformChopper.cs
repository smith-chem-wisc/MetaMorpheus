namespace EngineLayer.Truncation
{
    /// <summary>
    /// Pass 3 chopping algorithm: removes residues — together with any PTMs locked to them — one at a
    /// time from the indicated terminus of a parent proteoform until the precursor mass matches the
    /// scan within tolerance at an allowed notch. See 01_Architecture.md decisions #9, #10.
    ///
    /// Phase 0 stub: implemented in Phase 2.
    /// </summary>
    public static class ProteoformChopper
    {
        // TODO Phase 2: ChopUntilMassMatches(parent, terminusToChop, targetMass, tolerance, massDiffAcceptor)
        // returning the matched truncated form, or null when no integer-residue chop matches (#9).
    }
}
