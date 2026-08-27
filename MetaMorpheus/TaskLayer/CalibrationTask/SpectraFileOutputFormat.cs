namespace TaskLayer
{
    /// <summary>
    /// The container a task writes its processed spectra to.
    /// </summary>
    /// <remarks>
    /// mzML is the default and should stay that way. MGF carries a title, MS level, precursor m/z and
    /// charge, scan number, retention time and the peak list, and has nowhere to put dissociation type,
    /// analyzer type, scan filter, isolation window, precursor scan number or injection time. A task's
    /// output is fed to the next task in the chain, so choosing MGF means the search that follows sees a
    /// file without those fields.
    ///
    /// It is offered because "mgf in, mgf out" is a real workflow — smith-chem-wisc/MetaMorpheus#1673 —
    /// not because it is equivalent. Measured on the bundled yeast test file, a Calibrate then Search chain
    /// finds 85 PSMs through mzML and 61 through mgf.
    /// </remarks>
    public enum SpectraFileOutputFormat
    {
        MzML,
        Mgf
    }
}
