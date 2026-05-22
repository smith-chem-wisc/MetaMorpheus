namespace TaskLayer
{
    /// <summary>
    /// Parameters for <see cref="TruncationSearchTask"/>. Mirrors the SearchParameters pattern,
    /// holding only what the truncation search needs; expanded as later phases require more settings.
    /// Defaults encode the locked decisions in 01_Architecture.md (#3, #16, #17, #20).
    /// </summary>
    public class TruncationSearchParameters
    {
        /// <summary>
        /// Id of the upstream task whose deduped proteoform-level result seeds Pass 2 via the
        /// in-memory <see cref="TaskChainContext"/> hand-off (decision #1). Null for standalone runs.
        /// </summary>
        public string UpstreamSearchTaskId { get; set; } = null;

        /// <summary>
        /// Disk fallback (decision #1): path to an AllProteoforms.psmtsv to ingest when no upstream
        /// task result is present in the <see cref="TaskChainContext"/>.
        /// </summary>
        public string Pass1ProteoformsFilePath { get; set; } = null;

        /// <summary>
        /// Parent inclusion q-value threshold (decision #3). Deliberately permissive (10% FDR) to cast
        /// a wide net of parents; the truncation search's own FDR analysis polices the final output.
        /// </summary>
        public double ParentQValueThreshold { get; set; } = 0.10;

        /// <summary>
        /// Notch acceptor used when chopping a parent down to the precursor mass in Pass 3 (decision #80).
        /// Inherited from the Pass 1 search convention; defaults to the top-down ThreeMM notches so
        /// precursor isotope errors do not block a chop. Also sets the FDR notch count for the pooled
        /// analysis.
        /// </summary>
        public MassDiffAcceptorType MassDiffAcceptorType { get; set; } = MassDiffAcceptorType.ThreeMM;

        /// <summary>Custom mass-diff-acceptor string, used only when <see cref="MassDiffAcceptorType"/> is Custom.</summary>
        public string CustomMdac { get; set; } = null;

        /// <summary>Custom output-folder name used for benchmark runs (decision #20). Null = default naming.</summary>
        public string CustomOutputFolderName { get; set; } = null;

        /// <summary>Write decoy rows to the output TSVs (decision #17).</summary>
        public bool WriteDecoys { get; set; } = true;

        /// <summary>Write contaminant rows to the output TSVs (decision #17).</summary>
        public bool WriteContaminants { get; set; } = true;
    }
}
