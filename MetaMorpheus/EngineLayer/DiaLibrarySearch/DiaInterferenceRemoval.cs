#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using StatisticalModels;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>
/// Removal of interfering precursors, after DIA-NN (remove_ifs): a match whose most intense library fragments are already
/// explained by a better-scoring match co-eluting in the same or a neighbouring isolation window is dropped. Such a match
/// rides on the other's signal. Targets and decoys are treated alike, so a false target riding on a real peptide goes the
/// same way as a decoy would, and the target-decoy estimate stays fair.
/// </summary>
public static class DiaInterferenceRemoval
{
    /// <summary>A match's most intense library fragments that are checked against a better match.</summary>
    public const int CheckedFragments = 6;

    /// <summary>Only matches scoring at least as well as the target at this provisional q-value take part; the rest
    /// cannot be reported anyway, and leaving them out keeps the pass cheap.</summary>
    public const double ProvisionalQValue = 0.10;

    /// <param name="matches">One match per precursor, already scored.</param>
    /// <param name="fragmentMzsByIntensity">A precursor's library fragment m/z, most intense first.</param>
    /// <param name="windowOf">The index of the isolation window a precursor m/z falls in.</param>
    /// <param name="rtToleranceMinutes">How far apart two apexes may be and still co-elute.</param>
    /// <param name="explainedFragments">How many of the checked fragments a better match must explain.</param>
    /// <param name="sameMzOnly">
    /// As DIA-NN: only a better match with the same precursor m/z, or with it as its +1 isotope, can explain another.
    /// </param>
    /// <param name="ownIonCoElution">The co-elution of a winner's own ions against a rival it beat (<see cref="OwnIonCoElution"/>); called as (winner, rival).</param>
    /// <param name="ownIonFloor">When above 0, a winner below this against any rival it beat goes too. A value that is not a number is no evidence against it.</param>
    /// <returns>The matches that are not explained by a better one, in their input order.</returns>
    public static List<DiaPrecursorMatch> Remove(IReadOnlyList<DiaPrecursorMatch> matches, Func<int, float[]> fragmentMzsByIntensity,
        Func<double, int> windowOf, double rtToleranceMinutes, double tolerancePpm, int explainedFragments, Action<string>? report = null,
        bool sameMzOnly = false, Func<DiaPrecursorMatch, DiaPrecursorMatch, double>? ownIonCoElution = null, double ownIonFloor = 0)
    {
        ArgumentNullException.ThrowIfNull(matches);
        ArgumentNullException.ThrowIfNull(fragmentMzsByIntensity);
        ArgumentNullException.ThrowIfNull(windowOf);
        if (matches.Count == 0)
            return [];

        double[] q = TargetDecoyQValues.Compute(matches.Select(m => m.Score).ToArray(), matches.Select(m => m.IsDecoy).ToArray());
        double floor = matches.Where((m, i) => q[i] <= ProvisionalQValue).Select(m => m.Score).DefaultIfEmpty(double.PositiveInfinity).Min();

        var clock = System.Diagnostics.Stopwatch.StartNew();
        long neighbours = 0;
        int considered = 0;
        var removed = new HashSet<int>();
        var kept = new Dictionary<int, List<(double Rt, float[] Mzs, double PrecursorMz, int Charge, DiaPrecursorMatch Match)>>();
        bool askOwnIons = ownIonCoElution is not null && ownIonFloor > 0;
        var rivalsOf = new Dictionary<int, (DiaPrecursorMatch Winner, List<DiaPrecursorMatch> Rivals)>();
        foreach (var match in matches.Where((m, i) => q[i] <= ProvisionalQValue).OrderByDescending(m => m.Score).ThenBy(m => m.PrecursorIndex))
        {
            considered++;
            int window = windowOf(match.PrecursorMz);
            float[] own = fragmentMzsByIntensity(match.PrecursorIndex);
            float[] checkedMzs = own.Take(CheckedFragments).ToArray();

            // How many checked fragments a better match explains, or -1 when it cannot explain this one at all
            int ExplainedBy(double rt, float[] mzs, double betterMz, int betterCharge) =>
                Math.Abs(rt - match.ApexRt.Value) > rtToleranceMinutes || (sameMzOnly && !SameOrIsotope(match.PrecursorMz, betterMz, betterCharge, tolerancePpm))
                    ? -1
                    : checkedMzs.Count(mz => mzs.Any(other => Math.Abs(other - mz) <= mz * tolerancePpm * 1e-6));

            bool explained = false;
            for (int w = window - 1; w <= window + 1 && !explained; w++)
            {
                if (!kept.TryGetValue(w, out var better))
                    continue;
                foreach (var (rt, mzs, betterMz, betterCharge, _) in better)
                {
                    int count = ExplainedBy(rt, mzs, betterMz, betterCharge);
                    if (count < 0)
                        continue;
                    neighbours++;
                    if (count >= explainedFragments)
                    {
                        explained = true;
                        break;
                    }
                }
            }
            if (explained)
            {
                removed.Add(match.PrecursorIndex);
                if (askOwnIons)
                {
                    // Several better matches can explain it. Its rival is the one explaining the most, and the better of those on a tie
                    DiaPrecursorMatch? winner = null;
                    int most = 0;
                    for (int w = window - 1; w <= window + 1; w++)
                    {
                        if (!kept.TryGetValue(w, out var better))
                            continue;
                        foreach (var (rt, mzs, betterMz, betterCharge, candidate) in better)
                        {
                            int count = ExplainedBy(rt, mzs, betterMz, betterCharge);
                            if (count >= explainedFragments && (count > most || (count == most && candidate.Score > winner!.Score)))
                            {
                                most = count;
                                winner = candidate;
                            }
                        }
                    }
                    if (!rivalsOf.TryGetValue(winner!.PrecursorIndex, out var entry))
                        rivalsOf[winner.PrecursorIndex] = entry = (winner, []);
                    entry.Rivals.Add(match);
                }
                continue;
            }
            if (!kept.TryGetValue(window, out var list))
                kept[window] = list = [];
            list.Add((match.ApexRt.Value, own, match.PrecursorMz, match.Charge, match));
        }
        report?.Invoke($"interference removal: {considered} of {matches.Count} matches above the floor {floor:G4}; RT tolerance {rtToleranceMinutes:F4} min; " +
            $"{neighbours} co-eluting pairs checked; {removed.Count} removed ({matches.Count(m => m.IsDecoy && removed.Contains(m.PrecursorIndex))} decoys)  [{clock.Elapsed:mm\\:ss}]");
        if (askOwnIons)
        {
            // The rule keeps the rival removed: telling the two apart by these ions alone was worse than the score. A winner must stand against every rival it beat
            var dropped = rivalsOf.Values.Where(w => w.Rivals.Any(rival => ownIonCoElution!(w.Winner, rival) < ownIonFloor)).Select(w => w.Winner).ToList();
            removed.UnionWith(dropped.Select(m => m.PrecursorIndex));
            report?.Invoke($"interference removal: {dropped.Count} winners dropped, own ions below {ownIonFloor:G4} ({dropped.Count(m => m.IsDecoy)} decoys)");
        }
        return matches.Where(m => !removed.Contains(m.PrecursorIndex)).ToList();
    }

    /// <summary>
    /// How well a winner's own ions co-elute with what it shares with its rival: the library-weighted correlation of each own ion's trace with the summed trace of the shared
    /// ions. An ion is shared when the other precursor's whole library list holds that m/z, scorable or not, since that signal cannot tell the two apart. An own ion with no
    /// signal scores 0, so a strong library ion that is missing counts against the winner.
    /// </summary>
    /// <param name="winnerIons">The winner's scored fragments (m/z, library intensity).</param>
    /// <param name="winnerAllMzs">Every fragment m/z on the winner's library list.</param>
    /// <param name="rivalIons">The rival's scored fragments.</param>
    /// <param name="rivalAllMzs">Every fragment m/z on the rival's library list.</param>
    /// <param name="windows">The scans of each isolation window holding the precursor, each in retention-time order. The one with the most shared signal is measured.</param>
    /// <param name="rtMinutes">Where the two co-elute. The shared ions' apex is sought within twice the apex span of it.</param>
    /// <param name="tolerancePpm">The fragment tolerance.</param>
    /// <param name="apexHalfWidthScans">The half-width of the correlation around the shared apex.</param>
    /// <param name="mostIntensePeak">Take the most intense peak in tolerance rather than the closest, as the search does.</param>
    /// <returns>The correlation, 0 to 1; 0 when the winner has no own ions; not a number when the shared ions have no signal, so there is nothing to measure against.</returns>
    public static double OwnIonCoElution(IReadOnlyList<(double Mz, double Intensity)> winnerIons, IReadOnlyList<double> winnerAllMzs,
        IReadOnlyList<(double Mz, double Intensity)> rivalIons, IReadOnlyList<double> rivalAllMzs, IReadOnlyList<MsDataScan[]> windows,
        double rtMinutes, double tolerancePpm, int apexHalfWidthScans, bool mostIntensePeak = true)
    {
        bool Near(double a, double b) => Math.Abs(a - b) <= a * tolerancePpm * 1e-6;
        var own = winnerIons.Where(f => !rivalAllMzs.Any(o => Near(f.Mz, o))).ToList();
        var shared = winnerIons.Where(f => rivalAllMzs.Any(o => Near(f.Mz, o)))
            .Concat(rivalIons.Where(f => winnerAllMzs.Any(o => Near(f.Mz, o)))).Select(f => f.Mz).ToList();
        var tolerance = new PpmTolerance(tolerancePpm);
        int halfSpan = 2 * (apexHalfWidthScans + 1);

        double[] Trace(MsDataScan[] scans, int first, int count, double mz)
        {
            var trace = new double[count];
            for (int k = 0; k < count; k++)
            {
                var spectrum = scans[first + k].MassSpectrum;
                if (spectrum.Size == 0)
                    continue;
                int peak = DiaLibrarySearchEngine.FragmentPeakIndex(spectrum, mz, tolerance, mostIntensePeak);
                if (peak >= 0)
                    trace[k] = spectrum.YArray[peak];
            }
            return trace;
        }

        double bestSignal = 0;
        int bestApex = -1, firstScan = 0;
        MsDataScan[] bestWindow = [];
        double[] sharedTrace = [];
        foreach (var scans in windows)
        {
            if (scans.Length == 0)
                continue;
            int centre = Math.Min(FirstAtOrAfter(scans, rtMinutes), scans.Length - 1);
            int first = Math.Max(0, centre - halfSpan), last = Math.Min(scans.Length - 1, centre + halfSpan);
            var sum = new double[last - first + 1];
            foreach (double mz in shared)
            {
                var trace = Trace(scans, first, sum.Length, mz);
                for (int k = 0; k < sum.Length; k++)
                    sum[k] += trace[k];
            }
            int apex = Array.IndexOf(sum, sum.Max());
            if (sum[apex] > bestSignal)
            {
                bestSignal = sum[apex];
                bestApex = apex;
                bestWindow = scans;
                firstScan = first;
                sharedTrace = sum;
            }
        }
        if (bestApex < 0)
            return double.NaN;
        double total = own.Sum(f => f.Intensity);
        if (total <= 0)
            return 0;
        int from = Math.Max(0, bestApex - apexHalfWidthScans), to = Math.Min(sharedTrace.Length - 1, bestApex + apexHalfWidthScans);
        return own.Sum(f => f.Intensity * DiaLibrarySearchEngine.ClippedPearson(Trace(bestWindow, firstScan, sharedTrace.Length, f.Mz), sharedTrace, from, to)) / total;
    }

    /// <summary>The share of a winner's library intensity on ions its rival's whole list does not hold; 0 when it has no ions.</summary>
    public static double OwnIonShare(IReadOnlyList<(double Mz, double Intensity)> winnerIons, IReadOnlyList<double> rivalAllMzs, double tolerancePpm)
    {
        double total = winnerIons.Sum(f => f.Intensity);
        return total > 0
            ? winnerIons.Where(f => !rivalAllMzs.Any(o => Math.Abs(f.Mz - o) <= f.Mz * tolerancePpm * 1e-6)).Sum(f => f.Intensity) / total
            : 0;
    }

    /// <summary>The first scan at or after a retention time, or the scan count when every scan is earlier.</summary>
    private static int FirstAtOrAfter(MsDataScan[] scans, double rtMinutes)
    {
        int low = 0, high = scans.Length;
        while (low < high)
        {
            int middle = (low + high) / 2;
            if (scans[middle].RetentionTime < rtMinutes) low = middle + 1; else high = middle;
        }
        return low;
    }

    /// <summary>The weaker precursor's m/z is the better one's, or the better one's +1 isotope.</summary>
    private static bool SameOrIsotope(double mz, double betterMz, int betterCharge, double tolerancePpm) =>
        Math.Abs(mz - betterMz) <= mz * tolerancePpm * 1e-6
        || Math.Abs(mz - (betterMz + 1.0033548 / Math.Max(1, betterCharge))) <= mz * tolerancePpm * 1e-6;
}
