namespace EngineLayer.Truncation
{
    /// <summary>
    /// Pass 2 indexed scoring engine for the truncation search: builds a fragment index from the
    /// Pass 1 parent proteoforms and computes, per MS2 scan, two single-series scores
    /// (N-terminal-ion series only, C-terminal-ion series only) against every candidate parent,
    /// yielding the winning (parent, terminus) pair. See 01_Architecture.md decisions #4-#8.
    ///
    /// Phase 0 stub: implemented in Phase 1.
    /// </summary>
    public class TruncationSearchEngine
    {
        // TODO Phase 1: build the fragment index (following CommonParameters.DissociationType) and run
        // the per-scan dual single-series scoring loop, honoring the intact-skip rule (#4a).
    }
}
