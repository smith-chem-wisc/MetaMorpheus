#nullable enable

namespace EngineLayer.ParallelSearch.FdrAlignment;

public readonly record struct ProteinGroupBaselineFdrEntry : IScoreBaselineEntry
{
    public ProteinGroupBaselineFdrEntry(
        double score,
        double qValue,
        double bestPeptideQValue,
        double bestPeptidePep,
        int cumulativeTarget,
        int cumulativeDecoy)
    {
        Score = score;
        QValue = qValue;
        BestPeptideQValue = bestPeptideQValue;
        BestPeptidePep = bestPeptidePep;
        CumulativeTarget = cumulativeTarget;
        CumulativeDecoy = cumulativeDecoy;
    }

    public double Score { get; }

    public double QValue { get; }

    public double BestPeptideQValue { get; }

    public double BestPeptidePep { get; }

    public int CumulativeTarget { get; }

    public int CumulativeDecoy { get; }
}
