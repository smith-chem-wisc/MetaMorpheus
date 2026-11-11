#nullable enable

using System.Collections.Generic;
using System.Linq;
using EngineLayer.FdrAnalysis;
using EngineLayer.SpectrumMatch;

namespace EngineLayer.ParallelSearch.FdrAlignment;

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

/// <summary>
/// Gives transient PSMs, peptides and protein groups the FDR values of baseline items that rank the same,
/// so a transient database's results are judged against the baseline search's FDR instead of a new
/// calculation on a small subset. Each item is looked up on its own, so the order items arrive in does
/// not matter.
/// </summary>
public static class BaselineFdr
{
    public static BaselineScoreLookup<FdrInfo> CreatePsmLookup(IEnumerable<SpectralMatch> baselinePsms)
    {
        return new BaselineScoreLookup<FdrInfo>(baselinePsms
            .Where(psm => psm is not null)
            .Select(psm => (psm.Score, (psm.PsmFdrInfo ?? new FdrInfo()).Clone())));
    }

    /// <summary>
    /// Only baseline PSMs that were given peptide-level FDR take part.
    /// </summary>
    public static BaselineScoreLookup<FdrInfo> CreatePeptideLookup(IEnumerable<SpectralMatch> baselinePsms)
    {
        return new BaselineScoreLookup<FdrInfo>(baselinePsms
            .Where(psm => psm?.PeptideFdrInfo is not null && psm.PeptideFdrInfo.QValue < 2)
            .Select(psm => (psm.Score, psm.PeptideFdrInfo.Clone())));
    }

    /// <param name="filterType">The filter the baseline's protein FDR used (<see cref="ProteinScoringAndFdrEngine.FilterType"/>).</param>
    public static ProteinGroupFdrLookup CreateProteinGroupLookup(IEnumerable<ProteinGroup> baselineGroups, FilterType filterType)
    {
        return new ProteinGroupFdrLookup(baselineGroups
            .Where(group => group is not null)
            .Select(group => (ProteinGroupRank.Of(group, filterType),
                new ProteinGroupFdr(group.QValue, group.CumulativeTarget, group.CumulativeDecoy))),
            filterType);
    }

    public static void AlignPsms(BaselineScoreLookup<FdrInfo> lookup, IEnumerable<SpectralMatch> psms)
    {
        if (lookup.IsEmpty)
            return;

        foreach (var psm in psms)
            psm.PsmFdrInfo = lookup.Find(psm.Score).Clone();
    }

    public static void AlignPeptides(BaselineScoreLookup<FdrInfo> lookup, IEnumerable<SpectralMatch> peptides)
    {
        if (lookup.IsEmpty)
            return;

        foreach (var peptide in peptides)
            peptide.PeptideFdrInfo = lookup.Find(peptide.Score).Clone();
    }

    /// <summary>
    /// The groups' best-peptide values must already be set from their own PSMs
    /// (<see cref="ProteinScoringAndFdrEngine"/>'s PopulateBestPeptideMetrics); they decide where each group ranks.
    /// </summary>
    public static void AlignProteinGroups(ProteinGroupFdrLookup lookup, IEnumerable<ProteinGroup> groups)
    {
        if (lookup.IsEmpty)
            return;

        foreach (var group in groups)
        {
            var fdr = lookup.Find(ProteinGroupRank.Of(group, lookup.FilterType));
            group.QValue = fdr.QValue;
            group.CumulativeTarget = fdr.CumulativeTarget;
            group.CumulativeDecoy = fdr.CumulativeDecoy;
        }
    }
}
