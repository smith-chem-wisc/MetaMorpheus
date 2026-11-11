#nullable enable

using EngineLayer.FdrAnalysis;

namespace EngineLayer.ParallelSearch.FdrAlignment;

public sealed class PsmSpectralMatchFdrAlignmentService : ScoreBasedFdrAlignmentServiceBase<SpectralMatch, SpectralMatchBaselineFdrEntry>
{
    protected override double GetScore(SpectralMatch item)
    {
        return item.Score;
    }

    protected override bool TryBuildBaselineEntry(SpectralMatch item, out SpectralMatchBaselineFdrEntry entry)
    {
        entry = new SpectralMatchBaselineFdrEntry(item.Score, item.PsmFdrInfo ?? new FdrInfo());
        return true;
    }

    protected override void ApplyBaselineEntry(SpectralMatch item, SpectralMatchBaselineFdrEntry entry)
    {
        item.PsmFdrInfo = entry.FdrInfo.Clone();
    }
}
