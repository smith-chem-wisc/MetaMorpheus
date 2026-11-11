#nullable enable

using System;
using System.Collections.Generic;
using System.Linq;
using EngineLayer.FdrAnalysis;

namespace EngineLayer.ParallelSearch.FdrAlignment;

/// <summary>
/// Baseline values sorted by score, highest first. A transient score is given the value of the
/// lowest-scoring baseline entry that still scores at least as high; a score above every baseline
/// score gets the first entry, and one below every baseline score gets the last.
/// </summary>
public sealed class BaselineScoreLookup<T>
{
    private readonly double[] _scores;
    private readonly T[] _values;

    public BaselineScoreLookup(IEnumerable<(double Score, T Value)> baseline)
    {
        ArgumentNullException.ThrowIfNull(baseline);

        // OrderByDescending is stable, so equal scores keep the order they were given in
        var sorted = baseline.OrderByDescending(entry => entry.Score).ToArray();
        _scores = sorted.Select(entry => entry.Score).ToArray();
        _values = sorted.Select(entry => entry.Value).ToArray();
    }

    public int Count => _scores.Length;

    public bool IsEmpty => _scores.Length == 0;

    public T Find(double score)
    {
        if (IsEmpty)
            throw new InvalidOperationException("The baseline lookup is empty.");

        // count the baseline entries scoring at least as high as the transient score
        int low = 0;
        int high = _scores.Length;
        while (low < high)
        {
            int middle = (low + high) / 2;
            if (_scores[middle] >= score)
                low = middle + 1;
            else
                high = middle;
        }

        return _values[Math.Max(low - 1, 0)];
    }
}

/// <summary>
/// The FDR values a transient protein group takes from the baseline.
/// </summary>
public readonly record struct ProteinGroupFdr(
    double QValue,
    double BestPeptideQValue,
    double BestPeptidePep,
    int CumulativeTarget,
    int CumulativeDecoy);

/// <summary>
/// Gives transient PSMs, peptides and protein groups the FDR values of baseline items with the same
/// score, so a transient database's results are judged against the baseline search's FDR instead of
/// a new calculation. Each item is looked up on its own, so the order items arrive in does not matter.
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

    public static BaselineScoreLookup<ProteinGroupFdr> CreateProteinGroupLookup(IEnumerable<ProteinGroup> baselineGroups)
    {
        return new BaselineScoreLookup<ProteinGroupFdr>(baselineGroups
            .Where(group => group is not null)
            .Select(group => (group.ProteinGroupScore, new ProteinGroupFdr(
                group.QValue,
                group.BestPeptideQValue,
                group.BestPeptidePEP,
                group.CumulativeTarget,
                group.CumulativeDecoy))));
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

    public static void AlignProteinGroups(BaselineScoreLookup<ProteinGroupFdr> lookup, IEnumerable<ProteinGroup> groups)
    {
        if (lookup.IsEmpty)
            return;

        foreach (var group in groups)
        {
            var fdr = lookup.Find(group.ProteinGroupScore);
            group.QValue = fdr.QValue;
            group.BestPeptideQValue = fdr.BestPeptideQValue;
            group.BestPeptidePEP = fdr.BestPeptidePep;
            group.CumulativeTarget = fdr.CumulativeTarget;
            group.CumulativeDecoy = fdr.CumulativeDecoy;
        }
    }
}
