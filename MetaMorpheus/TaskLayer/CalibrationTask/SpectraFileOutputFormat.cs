namespace TaskLayer
{
    /// <summary>
    /// The container a task writes its processed spectra to.
    /// </summary>
    /// <remarks>
    /// mzML is the default and should stay that way. A task's output is fed to the next task in the chain,
    /// so this choice decides what the search that follows can see.
    ///
    /// Offered because "mgf in, mgf out" is a real workflow — smith-chem-wisc/MetaMorpheus#1673. Measured on
    /// the bundled yeast file, a Calibrate then Search chain produces the same 80 PSMs through Mgf as
    /// through MzML, and 63 through MgfMs2Only.
    /// </remarks>
    public enum SpectraFileOutputFormat
    {
        /// <summary>Default. Carries every field a later task reads.</summary>
        MzML,

        /// <summary>
        /// MGF including MS1 scans, so a later task can still deconvolute precursors. Out of spec: Matrix
        /// Science requires a PEPMASS in every block and an MS1 has none. ProteoWizard and OpenMS read it,
        /// MSToolkit — and so Comet — does not. Use <see cref="MgfMs2Only"/> to feed one of those.
        /// </summary>
        Mgf,

        /// <summary>
        /// MGF with MS2 and above only, which every reader accepts. No MS1 means no precursor
        /// deconvolution downstream, costing roughly a fifth of the identifications.
        /// </summary>
        MgfMs2Only
    }
}
