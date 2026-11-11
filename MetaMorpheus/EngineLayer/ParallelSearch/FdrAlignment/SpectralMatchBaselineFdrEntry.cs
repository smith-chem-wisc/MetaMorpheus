#nullable enable

using EngineLayer.FdrAnalysis;

namespace EngineLayer.ParallelSearch.FdrAlignment;

public readonly record struct SpectralMatchBaselineFdrEntry : IScoreBaselineEntry
{
    public double Score { get; }

    public FdrInfo FdrInfo { get; }

    public SpectralMatchBaselineFdrEntry(double score, FdrInfo psmFdrInfo)
    {
        Score = score;
        FdrInfo = psmFdrInfo.Clone();
    }
}
