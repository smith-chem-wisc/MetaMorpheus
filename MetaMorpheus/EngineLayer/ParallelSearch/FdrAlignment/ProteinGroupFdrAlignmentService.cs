#nullable enable

namespace EngineLayer.ParallelSearch.FdrAlignment;

public sealed class ProteinGroupFdrAlignmentService : ScoreBasedFdrAlignmentServiceBase<ProteinGroup, ProteinGroupBaselineFdrEntry>
{
    protected override double GetScore(ProteinGroup item)
    {
        return item.ProteinGroupScore;
    }

    protected override bool TryBuildBaselineEntry(ProteinGroup item, out ProteinGroupBaselineFdrEntry entry)
    {
        entry = new ProteinGroupBaselineFdrEntry(
            item.ProteinGroupScore,
            item.QValue,
            item.BestPeptideQValue,
            item.BestPeptidePEP,
            item.CumulativeTarget,
            item.CumulativeDecoy);
        return true;
    }

    protected override void ApplyBaselineEntry(ProteinGroup item, ProteinGroupBaselineFdrEntry entry)
    {
        item.QValue = entry.QValue;
        item.BestPeptideQValue = entry.BestPeptideQValue;
        item.BestPeptidePEP = entry.BestPeptidePep;
        item.CumulativeTarget = entry.CumulativeTarget;
        item.CumulativeDecoy = entry.CumulativeDecoy;
    }
}
