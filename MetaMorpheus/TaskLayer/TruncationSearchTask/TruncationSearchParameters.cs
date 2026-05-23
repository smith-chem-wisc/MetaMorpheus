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

        /// <summary>
        /// Optional path to the rolling, append-only <c>perf_log.tsv</c> (03_Benchmarks.md). When set, the
        /// task appends one metrics+timing row per run. Null (default) = no perf logging, so ordinary runs
        /// and tests are unaffected.
        /// </summary>
        public string PerfLogPath { get; set; } = null;

        /// <summary>
        /// Seed truncation parents from the DATABASE — theoretical full-length proteoforms generated like a
        /// top-down search (digest each protein with the configured fixed/variable mods, so annotated PTMs
        /// ride along) — rather than only from upstream-identified proteoforms. Removes the requirement that a
        /// truncation's parent be observed intact, at the cost of a much larger Pass 2 search. Single-terminus
        /// only (internal truncations are not generated here). Off by default.
        /// </summary>
        public bool SeedParentsFromDatabase { get; set; } = false;

        /// <summary>
        /// When database-seeding (<see cref="SeedParentsFromDatabase"/>), first narrow the database to
        /// proteins supported by de-novo sequence tags extracted from the scans (doc §11.2.4) rather than
        /// digesting every protein. Global-union variant: the candidate set is the union over all scans, so
        /// the indexed engine is reused unchanged — this controls the parent explosion (the all-DB blowup)
        /// but does not yet restrict candidates per scan. Off by default.
        /// </summary>
        public bool UseSequenceTagFilter { get; set; } = false;

        /// <summary>
        /// With the sequence-tag filter on, score each scan ONLY against parents from its own tag-selected
        /// proteins (per-scan restriction) rather than the global union. Gives per-scan selectivity and stops a
        /// scan being won by a protein its own tags did not support (the winner-competition loss). Off = global
        /// union. Requires <see cref="UseSequenceTagFilter"/>.
        /// </summary>
        public bool UsePerScanTagRestriction { get; set; } = false;

        /// <summary>Length (residues) of the de-novo sequence tags and the matching protein k-mer index.</summary>
        public int TagLength { get; set; } = 4;

        /// <summary>Minimum number of distinct extracted tags a protein must contain to be a tag-filter candidate.</summary>
        public int MinTagHits { get; set; } = 2;

        /// <summary>
        /// Maximum parent proteoform monoisotopic mass (Da) placed in the Pass 2 index (decision #6); heavier
        /// parents are excluded. Defaults to the search convention (30 kDa). RAISE this for database-seeded
        /// runs so large proteins — whose small fragments populate an LMW fraction but whose intact form is
        /// never observed — can still be chopped down to a matching truncation.
        /// </summary>
        public double MaxParentMass { get; set; } = EngineLayer.Truncation.TruncationSearchEngine.DefaultMaxFragmentSize;

        /// <summary>Write decoy rows to the output TSVs (decision #17).</summary>
        public bool WriteDecoys { get; set; } = true;

        /// <summary>Write contaminant rows to the output TSVs (decision #17).</summary>
        public bool WriteContaminants { get; set; } = true;
    }
}
