#nullable enable

using System.Collections.Generic;
using EngineLayer.SpectrumMatch;

namespace EngineLayer.ParallelSearch.FdrAlignment;

public sealed class ProteinGroupFdrLookup : BaselineRankLookup<ProteinGroupRank, ProteinGroupFdr>
{
    public ProteinGroupFdrLookup(IEnumerable<(ProteinGroupRank Rank, ProteinGroupFdr Fdr)> baseline, FilterType filterType)
        : base(baseline, ProteinGroupRank.BetterFirst)
    {
        FilterType = filterType;
    }

    /// <summary>
    /// The filter the baseline's protein FDR ranked groups by. Transient groups must be ranked the same way.
    /// </summary>
    public FilterType FilterType { get; }
}

/// <summary>
/// Where a protein group sits in protein FDR's ranking: its best peptide's q-value (or PEP, when FDR filters
/// on PEP q-value), lower first, then its best peptide's score, higher first. This is the order
/// <see cref="ProteinScoringAndFdrEngine"/> assigns protein q-values in; the summed protein-group score is not.
/// </summary>
public readonly record struct ProteinGroupRank(double BestPeptideMetric, double BestPeptideScore)
{
    public static ProteinGroupRank Of(ProteinGroup group, FilterType filterType)
    {
        return new ProteinGroupRank(
            filterType == FilterType.PepQValue ? group.BestPeptidePEP : group.BestPeptideQValue,
            group.BestPeptideScore);
    }

    public static IComparer<ProteinGroupRank> BetterFirst { get; } = Comparer<ProteinGroupRank>.Create((a, b) =>
    {
        int byMetric = a.BestPeptideMetric.CompareTo(b.BestPeptideMetric);
        return byMetric != 0 ? byMetric : b.BestPeptideScore.CompareTo(a.BestPeptideScore);
    });
}

/// <summary>
/// The FDR values a transient protein group takes from the baseline. Its best-peptide values are its own.
/// </summary>
public readonly record struct ProteinGroupFdr(double QValue, int CumulativeTarget, int CumulativeDecoy);
