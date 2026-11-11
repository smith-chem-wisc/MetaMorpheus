#nullable enable

namespace EngineLayer.ParallelSearch.FdrAlignment;

public sealed class PeptideSpectralMatchFdrAlignmentService : ScoreBasedFdrAlignmentServiceBase<SpectralMatch, SpectralMatchBaselineFdrEntry>
{
    protected override double GetScore(SpectralMatch item)
    {
        return item.Score;
    }

    protected override bool TryBuildBaselineEntry(SpectralMatch item, out SpectralMatchBaselineFdrEntry entry)
    {
        if (item.PeptideFdrInfo is not null && item.PeptideFdrInfo.QValue < 2)
        {
            entry = new SpectralMatchBaselineFdrEntry(item.Score, item.PeptideFdrInfo);
            return true;
        }

        entry = default;
        return false;
    }

    protected override void ApplyBaselineEntry(SpectralMatch item, SpectralMatchBaselineFdrEntry entry)
    {
        item.PeptideFdrInfo = entry.FdrInfo.Clone();
    }
}
