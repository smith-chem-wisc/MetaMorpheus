namespace EngineLayer.Truncation
{
    /// <summary>
    /// Mutable accumulator for Pass 3 wall-clock timings (perf logging, docs/Truncation-Search.md, Benchmarking hook). Passed into
    /// <see cref="TruncationPass3.ScoreTruncation"/> across all winners so the chopping and standard-scoring
    /// portions are tallied separately. Null = timing not collected (default for normal runs and tests).
    /// </summary>
    public class TruncationTimings
    {
        public double ChoppingSeconds { get; set; }
        public double ScoringSeconds { get; set; }
    }
}
