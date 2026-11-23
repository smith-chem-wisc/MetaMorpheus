#nullable enable

using System.Collections.Generic;
using System.Linq;
using EngineLayer.FdrAnalysis;
using EngineLayer.SpectrumMatch;

namespace EngineLayer.ParallelSearch.FdrAlignment;

/// <summary>The notch q-value of one precursor notch, with the cumulative counts it came from.</summary>
public readonly record struct NotchFdr(double QValue, double CumulativeTarget, double CumulativeDecoy);

/// <summary>
/// Gives transient PSMs, peptides and protein groups the FDR values of baseline items that rank the same,
/// so a transient database's results are judged against the baseline search's FDR instead of a new
/// calculation on a small subset. Each item is looked up on its own, so the order items arrive in does
/// not matter.
/// <para>
/// With more than one precursor notch (e.g. OneMM) a match also has a notch q-value, which decides whether it counts as
/// confident (<see cref="FilteredPsms"/>, ParallelSearch's ConfidenceCutoff). It follows its own notch's curve, not the
/// score curve the q-value follows, so it is looked up per notch (<see cref="CreateNotchLookups"/>). Only match-level
/// values are written: a transient match shares its baseline hypotheses with the baseline match it was copied from.
/// </para>
/// </summary>
public static class BaselineFdr
{
    public static BaselineScoreLookup<FdrInfo> CreatePsmLookup(IEnumerable<SpectralMatch> baselinePsms)
    {
        return new BaselineScoreLookup<FdrInfo>(WithMonotonePepQValues(baselinePsms
            .Where(psm => psm is not null)
            .Select(psm => (psm.Score, (psm.PsmFdrInfo ?? new FdrInfo()).Clone()))));
    }

    /// <summary>
    /// Only baseline PSMs that were given peptide-level FDR take part.
    /// </summary>
    public static BaselineScoreLookup<FdrInfo> CreatePeptideLookup(IEnumerable<SpectralMatch> baselinePsms)
    {
        return new BaselineScoreLookup<FdrInfo>(WithMonotonePepQValues(baselinePsms
            .Where(psm => psm?.PeptideFdrInfo is not null && psm.PeptideFdrInfo.QValue < 2)
            .Select(psm => (psm.Score, psm.PeptideFdrInfo.Clone()))));
    }

    /// <summary>
    /// Replaces the PEP q-values of a score-ordered lookup with their isotonic fit: non-decreasing as the score falls
    /// (pool adjacent violators, ties in score pooled first).
    /// <para>
    /// A transient match takes the FDR values of the baseline entry at its score. The q-value is monotone in score, but the
    /// PEP q-value is not (it ranks by PEP, which uses more than the score), so under PEP filtering a transient match took
    /// whatever PEP q-value one neighbouring baseline PSM happened to have. The fit gives it the PEP q-value the
    /// baseline shows at that score. Unset PEP q-values (2) are left out and left alone.
    /// </para>
    /// </summary>
    internal static List<(double Score, FdrInfo Value)> WithMonotonePepQValues(IEnumerable<(double Score, FdrInfo Value)> entries)
    {
        var list = entries.ToList();
        var set = list.Where(e => e.Value.PEP_QValue < 2).OrderByDescending(e => e.Score).ToList();
        if (set.Count < 2)
            return list;

        // blocks of (sum, count, members); members in score order, one block per distinct score to start
        var sums = new List<double>();
        var counts = new List<int>();
        var members = new List<List<FdrInfo>>();
        double? lastScore = null;
        foreach (var (score, info) in set)
        {
            if (lastScore is double s && s == score)
            {
                sums[^1] += info.PEP_QValue;
                counts[^1]++;
                members[^1].Add(info);
            }
            else
            {
                sums.Add(info.PEP_QValue);
                counts.Add(1);
                members.Add([info]);
            }
            lastScore = score;

            // merge while the block before has a larger mean than this one (the fit must not fall as the score falls)
            while (sums.Count > 1 && sums[^2] / counts[^2] > sums[^1] / counts[^1])
            {
                sums[^2] += sums[^1];
                counts[^2] += counts[^1];
                members[^2].AddRange(members[^1]);
                sums.RemoveAt(sums.Count - 1);
                counts.RemoveAt(counts.Count - 1);
                members.RemoveAt(members.Count - 1);
            }
        }

        for (int b = 0; b < members.Count; b++)
        {
            double fitted = sums[b] / counts[b];
            foreach (var info in members[b])
                info.PEP_QValue = fitted;
        }

        return list;
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

    /// <summary>
    /// FDR gives a match whose tied hypotheses sit at different notches a notch q-value per hypothesis only, and leaves
    /// the match's own at its unset value of 2, which no cutoff passes. This gives each such baseline match the values of
    /// its best (lowest notch q-value) hypothesis, at the PSM level and, where peptide FDR ran, the peptide level: the
    /// values the PSM writer already reports for it. Run once after baseline FDR, before the lookups are made.
    /// </summary>
    public static void UseBestNotch(IEnumerable<SpectralMatch> baselineMatches)
    {
        foreach (var match in baselineMatches)
        {
            if (match is null || match.Notch is not null)
                continue;

            if (match.PsmFdrInfo is not null && BestHypothesisNotch(match, peptideLevel: false) is NotchFdr psmBest)
                SetNotch(match.PsmFdrInfo, psmBest);
            if (match.PeptideFdrInfo is not null && BestHypothesisNotch(match, peptideLevel: true) is NotchFdr peptideBest)
                SetNotch(match.PeptideFdrInfo, peptideBest);
        }
    }

    /// <summary>
    /// One score lookup per notch, from the baseline matches FDR reached at this level: a notch-unambiguous match
    /// contributes its own notch values, a notch-ambiguous one those of each notch its hypotheses sit at.
    /// </summary>
    public static IReadOnlyDictionary<int, BaselineScoreLookup<NotchFdr>> CreateNotchLookups(
        IEnumerable<SpectralMatch> baselineMatches, bool peptideLevel)
    {
        var entries = new Dictionary<int, List<(double Score, NotchFdr Value)>>();
        void Add(int notch, double score, NotchFdr value)
        {
            if (!entries.TryGetValue(notch, out var list))
                entries[notch] = list = new List<(double, NotchFdr)>();
            list.Add((score, value));
        }

        foreach (var match in baselineMatches)
        {
            var fdr = match?.GetFdrInfo(peptideLevel);
            if (fdr is null || fdr.QValue >= 2) // FDR did not reach it
                continue;

            if (match!.Notch is int notch)
            {
                Add(notch, match.Score, new NotchFdr(fdr.QValueNotch, fdr.CumulativeTargetNotch, fdr.CumulativeDecoyNotch));
                continue;
            }

            foreach (var group in match.BestMatchingBioPolymersWithSetMods.GroupBy(h => h.Notch))
            {
                if (HypothesisNotch(group.First(), peptideLevel) is NotchFdr value)
                    Add(group.Key, match.Score, value);
            }
        }

        return entries.ToDictionary(e => e.Key, e => new BaselineScoreLookup<NotchFdr>(e.Value));
    }

    /// <param name="notchLookups">With more than one notch, the baseline's per-notch lookups (<see cref="CreateNotchLookups"/>):
    /// each PSM then gets the best notch values at its score among the notches its hypotheses sit at.</param>
    public static void AlignPsms(BaselineScoreLookup<FdrInfo> lookup, IEnumerable<SpectralMatch> psms,
        IReadOnlyDictionary<int, BaselineScoreLookup<NotchFdr>>? notchLookups = null)
    {
        if (lookup.IsEmpty)
            return;

        foreach (var psm in psms)
        {
            psm.PsmFdrInfo = lookup.Find(psm.Score).Clone();
            AlignNotch(notchLookups, psm, psm.PsmFdrInfo);
        }
    }

    /// <param name="notchLookups">As for <see cref="AlignPsms"/>, at the peptide level.</param>
    public static void AlignPeptides(BaselineScoreLookup<FdrInfo> lookup, IEnumerable<SpectralMatch> peptides,
        IReadOnlyDictionary<int, BaselineScoreLookup<NotchFdr>>? notchLookups = null)
    {
        if (lookup.IsEmpty)
            return;

        foreach (var peptide in peptides)
        {
            peptide.PeptideFdrInfo = lookup.Find(peptide.Score).Clone();
            AlignNotch(notchLookups, peptide, peptide.PeptideFdrInfo);
        }
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

    /// <summary>A notch the baseline never reached leaves the values copied from the score lookup.</summary>
    private static void AlignNotch(IReadOnlyDictionary<int, BaselineScoreLookup<NotchFdr>>? notchLookups, SpectralMatch match,
        FdrInfo fdr)
    {
        if (notchLookups is null || notchLookups.Count == 0)
            return;

        NotchFdr? best = null;
        foreach (int notch in match.BestMatchingBioPolymersWithSetMods.Select(h => h.Notch).Distinct())
        {
            if (!notchLookups.TryGetValue(notch, out var notchLookup) || notchLookup.IsEmpty)
                continue;
            var value = notchLookup.Find(match.Score);
            if (best is null || value.QValue < best.Value.QValue)
                best = value;
        }

        if (best is NotchFdr found)
            SetNotch(fdr, found);
    }

    private static NotchFdr? BestHypothesisNotch(SpectralMatch match, bool peptideLevel)
    {
        NotchFdr? best = null;
        foreach (var hypothesis in match.BestMatchingBioPolymersWithSetMods)
        {
            if (HypothesisNotch(hypothesis, peptideLevel) is not NotchFdr value)
                return null; // FDR has not reached every hypothesis
            if (best is null || value.QValue < best.Value.QValue)
                best = value;
        }
        return best;
    }

    private static NotchFdr? HypothesisNotch(SpectralMatchHypothesis hypothesis, bool peptideLevel)
    {
        if (peptideLevel)
        {
            return hypothesis is { PeptideQValueNotch: double q, PeptideCumulativeTargetNotch: double target, PeptideCumulativeDecoyNotch: double decoy }
                ? new NotchFdr(q, target, decoy)
                : null;
        }

        return hypothesis is { QValueNotch: double psmQ, CumulativeTargetNotch: double psmTarget, CumulativeDecoyNotch: double psmDecoy }
            ? new NotchFdr(psmQ, psmTarget, psmDecoy)
            : null;
    }

    private static void SetNotch(FdrInfo fdr, NotchFdr value)
    {
        fdr.QValueNotch = value.QValue;
        fdr.CumulativeTargetNotch = value.CumulativeTarget;
        fdr.CumulativeDecoyNotch = value.CumulativeDecoy;
    }
}
