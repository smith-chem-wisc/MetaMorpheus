#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MassSpectrometry;
using MzLibUtil;
using StatisticalModels;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>
/// How much of a precursor's fragment signal confident precursors eluting with it already explain, and what is left of its evidence once that is taken out. A search scores each
/// precursor alone, so a false match is free to read signal that belongs to something else: decoys and false targets in a real target's score band share about twice as many
/// ions with what elutes beside them as they do with precursors far away in retention time, while targets share no more than chance. Once the classifier has given a first order,
/// that order says who the confident precursors are, and trees re-rank everything near the threshold on the score and these columns.
/// </summary>
/// <remarks>
/// Targets and decoys are measured alike against the same confident list, and only targets are on it, as a search knows nothing else. Each row is scored by trees that never saw its
/// group. As with <see cref="DiaTreeBlend"/>, the step decides only the order: the values stay the classifier's own, so later stages see the same score scale.
/// </remarks>
public static class DiaNeighbourEvidence
{
    /// <summary>Who counts as confident: targets at or under this provisional q-value.</summary>
    public const double ConfidentQValue = 0.01;

    /// <summary>Rows at or under this provisional q-value keep their places. Trees are coarse at the very top: unguarded, a few decoys climbing there cost 28% of the targets at q 0.1% on PXD005573 1 h.</summary>
    public const double KeptQValue = 0.001;

    /// <summary>The most a neighbour's own strongest ions are read to say how much of it is here.</summary>
    private const int ScalingIons = 6;

    /// <summary>With fewer of its top ions the row does not hold, a neighbour is scaled on all of them. A near-isomer has almost none of its own.</summary>
    private const int MinimumFreeIons = 3;

    /// <summary>An ion is clean while neighbours explain less than this share of it.</summary>
    private const double CleanBelow = 0.5;

    private const int StackFolds = 5;

    private static readonly DiaTreeBlend.Trees StackTrees = new(Count: 120, Leaves: 15, MinimumPerLeaf: 50, LearningRate: 0.05, FeatureFraction: 1, BaggingFraction: 1);

    /// <summary>The columns of <see cref="Measure"/>, in order.</summary>
    public static readonly string[] ColumnNames = ["ExplainedShare", "ExplainedShareTop6", "ResidualCosine", "CleanIons", "CleanLibraryShare", "CleanCoElution", "SharedIons", "Neighbours"];

    /// <summary>A confident precursor as a neighbour.</summary>
    /// <param name="ApexRt">Where it elutes, in minutes.</param>
    /// <param name="BaseSequence">Its sequence: a row's own charge states are not neighbours, since their signal is its own.</param>
    /// <param name="Mzs">Its scored library ions, most intense first.</param>
    /// <param name="Intensities">Their library intensities, in that order.</param>
    public sealed record Neighbour(double ApexRt, string BaseSequence, double[] Mzs, double[] Intensities);

    /// <summary>What the library says of a precursor.</summary>
    /// <param name="Ions">Its scored library ions, most intense first.</param>
    /// <param name="AllMzs">Every fragment m/z on its library list, scored or not.</param>
    /// <param name="BaseSequence">Its sequence.</param>
    public sealed record PrecursorIons(IReadOnlyList<(double Mz, double Intensity)> Ions, IReadOnlyList<double> AllMzs, string BaseSequence);

    /// <summary>
    /// Reads one row at its apex. Each neighbour eluting within the tolerance predicts its library intensity at each ion it shares with the row, scaled by the median of observed over
    /// library at its strongest ions the row does not hold. Three apex scans are summed. Always finite: a row with nothing to read has zeros.
    /// </summary>
    /// <param name="ions">The row's scored ions, most intense first.</param>
    /// <param name="allMzs">Every fragment m/z on the row's library list.</param>
    /// <param name="baseSequence">The row's sequence.</param>
    /// <param name="neighbours">The confident precursors of the row's isolation window, by apex retention time.</param>
    /// <param name="scans">The window's scans, in retention-time order.</param>
    /// <param name="apexScan">The row's apex among them.</param>
    /// <param name="rtToleranceMinutes">How far a neighbour's apex may be from the row's and still co-elute.</param>
    /// <param name="tolerancePpm">The fragment tolerance.</param>
    /// <param name="apexHalfWidthScans">The half-width of the co-elution around the apex.</param>
    /// <param name="mostIntensePeak">Take the most intense peak in tolerance rather than the closest, as the search does.</param>
    /// <returns>The columns named by <see cref="ColumnNames"/>.</returns>
    public static double[] Measure(IReadOnlyList<(double Mz, double Intensity)> ions, IReadOnlyList<double> allMzs, string baseSequence, IReadOnlyList<Neighbour> neighbours,
        MsDataScan[] scans, int apexScan, double rtToleranceMinutes, double tolerancePpm, int apexHalfWidthScans, bool mostIntensePeak = true)
    {
        var tolerance = new PpmTolerance(tolerancePpm);
        bool Near(double a, double b) => Math.Abs(a - b) <= a * tolerancePpm * 1e-6;
        double PeakAt(MsDataScan scan, double mz)
        {
            var spectrum = scan.MassSpectrum;
            if (spectrum.Size == 0)
                return 0;
            int peak = DiaLibrarySearchEngine.FragmentPeakIndex(spectrum, mz, tolerance, mostIntensePeak);
            return peak >= 0 ? spectrum.YArray[peak] : 0;
        }
        double AtApex(double mz)
        {
            double sum = 0;
            for (int k = Math.Max(0, apexScan - 1); k <= Math.Min(scans.Length - 1, apexScan + 1); k++)
                sum += PeakAt(scans[k], mz);
            return sum;
        }

        int count = ions.Count;
        var observed = new double[count];
        for (int i = 0; i < count; i++)
            observed[i] = AtApex(ions[i].Mz);
        double total = observed.Sum();

        double apexRt = scans[apexScan].RetentionTime;
        var predicted = new double[count];
        var shared = new bool[count];
        int coEluting = 0;
        for (int n = FirstAtOrAfter(neighbours, apexRt - rtToleranceMinutes); n < neighbours.Count && neighbours[n].ApexRt <= apexRt + rtToleranceMinutes; n++)
        {
            var neighbour = neighbours[n];
            if (neighbour.BaseSequence == baseSequence)
                continue;
            coEluting++;
            double[]? holds = null;
            for (int a = 0; a < count; a++)
                for (int b = 0; b < neighbour.Mzs.Length; b++)
                    if (Near(ions[a].Mz, neighbour.Mzs[b]))
                    {
                        (holds ??= new double[count])[a] += neighbour.Intensities[b];
                        shared[a] = true;
                    }
            if (holds is null)
                continue;

            int[] top = Enumerable.Range(0, Math.Min(ScalingIons, neighbour.Mzs.Length)).ToArray();
            int[] free = top.Where(j => !allMzs.Any(o => Near(o, neighbour.Mzs[j]))).ToArray();
            double[] ratios = (free.Length >= MinimumFreeIons ? free : top).Select(j => AtApex(neighbour.Mzs[j]) / Math.Max(1e-6, neighbour.Intensities[j])).Order().ToArray();
            double scale = ratios[ratios.Length / 2];
            for (int a = 0; a < count; a++)
                predicted[a] += scale * holds[a];
        }

        var taken = new double[count];
        var left = new double[count];
        for (int i = 0; i < count; i++)
        {
            taken[i] = Math.Min(observed[i], predicted[i]);
            left[i] = observed[i] - taken[i];
        }
        int topCount = Math.Min(6, count);
        double topTotal = observed.Take(topCount).Sum();
        double libraryTotal = ions.Sum(f => f.Intensity);
        int[] clean = Enumerable.Range(0, count).Where(i => observed[i] > 0 && taken[i] / observed[i] < CleanBelow).ToArray();

        return
        [
            total > 0 ? taken.Sum() / total : 0,
            topTotal > 0 ? taken.Take(topCount).Sum() / topTotal : 0,
            Cosine(ions, left),
            clean.Length,
            libraryTotal > 0 ? clean.Sum(i => ions[i].Intensity) / libraryTotal : 0,
            CoElution(clean),
            shared.Count(s => s),
            coEluting,
        ];

        // Each ion against the sum of the others, weighted by the library, so a lone ion has nothing to co-elute with and reads 0
        double CoElution(int[] set)
        {
            int first =Math.Max(0, apexScan - apexHalfWidthScans), last = Math.Min(scans.Length - 1, apexScan + apexHalfWidthScans), span = last - first + 1;
            var traces = new Dictionary<int, double[]>();
            var sum = new double[span];
            foreach (int i in set)
            {
                var trace = new double[span];
                for (int k = 0; k < span; k++)
                {
                    trace[k] = PeakAt(scans[first + k], ions[i].Mz);
                    sum[k] += trace[k];
                }
                traces[i] = trace;
            }
            double weighted = 0, weight = 0;
            var others = new double[span];
            foreach (int i in set)
            {
                for (int k = 0; k < span; k++)
                    others[k] = sum[k] - traces[i][k];
                double r = DiaLibrarySearchEngine.ClippedPearson(traces[i], others, 0, span - 1);
                weighted += ions[i].Intensity * (double.IsNaN(r) ? 0 : r);
                weight += ions[i].Intensity;
            }
            return weight > 0 ? weighted / weight : 0;
        }
    }

    /// <summary>Measures every row down to the provisional floor, each in the window holding its apex, against that window's confident targets.</summary>
    /// <param name="matches">One match per precursor, already scored.</param>
    /// <param name="ionsOf">What the library says of a precursor, by index; null when it cannot say.</param>
    /// <param name="windowBounds">Each isolation window's bounds.</param>
    /// <param name="windowScans">Each window's scans, in retention-time order.</param>
    /// <param name="rtToleranceMinutes">How far two apexes may be and still co-elute.</param>
    /// <param name="tolerancePpm">The fragment tolerance.</param>
    /// <param name="apexHalfWidthScans">The half-width of the co-elution around the apex.</param>
    /// <param name="mostIntensePeak">Take the most intense peak in tolerance rather than the closest, as the search does.</param>
    /// <param name="topFragments">How many of a row's ions are read. Neighbours hold all theirs.</param>
    /// <param name="threads">The most threads to read with.</param>
    /// <returns>Each match's columns; null for a match below the floor; zeros for one with no window or nothing on the library.</returns>
    public static double[]?[] MeasureAll(IReadOnlyList<DiaPrecursorMatch> matches, Func<int, PrecursorIons?> ionsOf, IReadOnlyList<(double Low, double High)> windowBounds,
        IReadOnlyList<MsDataScan[]> windowScans, double rtToleranceMinutes, double tolerancePpm, int apexHalfWidthScans, bool mostIntensePeak, int topFragments, int threads = 1)
    {
        var columns = new double[]?[matches.Count];
        if (matches.Count == 0)
            return columns;
        double[] q = TargetDecoyQValues.Compute(matches.Select(m => m.Score).ToArray(), matches.Select(m => m.IsDecoy).ToArray());
        double floor = matches.Where((_, i) => q[i] <= DiaInterferenceRemoval.ProvisionalQValue).Select(m => m.Score).DefaultIfEmpty(double.PositiveInfinity).Min();

        var confident = Enumerable.Range(0, matches.Count).Where(i => !matches[i].IsDecoy && q[i] <= ConfidentQValue).AsParallel().AsOrdered().WithDegreeOfParallelism(threads)
            .Select(i => (Match: matches[i], Ions: ionsOf(matches[i].PrecursorIndex))).Where(c => c.Ions is not null).ToList();
        var neighboursOf = windowBounds.Select(window => confident.Where(c => c.Match.PrecursorMz >= window.Low && c.Match.PrecursorMz < window.High)
            .OrderBy(c => c.Match.ApexRt.Value).ThenBy(c => c.Match.PrecursorIndex)
            .Select(c => new Neighbour(c.Match.ApexRt.Value, c.Ions!.BaseSequence, c.Ions.Ions.Select(f => f.Mz).ToArray(), c.Ions.Ions.Select(f => f.Intensity).ToArray())).ToArray()).ToArray();

        var options = new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, threads) };
        Parallel.For(0, matches.Count, options, i =>
        {
            var match = matches[i];
            if (!(match.Score >= floor))
                return;
            int own = -1, ownScan = -1;
            double gap = double.PositiveInfinity;
            for (int w = 0; w < windowBounds.Count; w++)
            {
                if (match.PrecursorMz < windowBounds[w].Low || match.PrecursorMz >= windowBounds[w].High || windowScans[w].Length == 0)
                    continue;
                int scan = NearestScan(windowScans[w], match.ApexRt.Value);
                double distance = Math.Abs(windowScans[w][scan].RetentionTime - match.ApexRt.Value);
                if (distance < gap)
                {
                    gap = distance;
                    own = w;
                    ownScan = scan;
                }
            }
            var library = own >= 0 ? ionsOf(match.PrecursorIndex) : null;
            columns[i] = library is null
                ? new double[ColumnNames.Length]
                : Measure(library.Ions.Take(topFragments).ToList(), library.AllMzs, library.BaseSequence, neighboursOf[own], windowScans[own], ownScan, rtToleranceMinutes, tolerancePpm, apexHalfWidthScans, mostIntensePeak);
        });
        return columns;
    }

    /// <summary>
    /// Re-orders rows by out-of-fold trees learning target against decoy from the score and the columns.
    /// </summary>
    /// <param name="groupKeys">Rows sharing a key, such as a peptide's charge states, stay in one fold.</param>
    /// <param name="keepAtOrAbove">Rows scoring at least this keep their places. The trees still learn from them.</param>
    /// <returns>The same score values as <paramref name="scores"/>, reordered, or a copy of it when there is too little to learn from.</returns>
    /// <exception cref="ArgumentException">The lists differ in length.</exception>
    public static double[] Stack(IReadOnlyList<double> scores, IReadOnlyList<double[]> columns, IReadOnlyList<bool> isDecoy, IReadOnlyList<string> groupKeys, int seed, int threads = 1,
        double keepAtOrAbove = double.PositiveInfinity)
    {
        int n = scores.Count;
        if (columns.Count != n || isDecoy.Count != n || groupKeys.Count != n)
            throw new ArgumentException("Scores, columns, decoy flags and group keys must have one entry per row.");
        double[][] features = scores.Select((score, i) => columns[i].Prepend(score).ToArray()).ToArray();
        double[]? margins = DiaTreeBlend.OutOfFoldMargins(features, scores, isDecoy, groupKeys, n, StackFolds, StackTrees, seed, threads);
        if (margins is null)
            return scores.ToArray();
        int[] order = Enumerable.Range(0, n).OrderByDescending(i => scores[i] >= keepAtOrAbove)
            .ThenByDescending(i => scores[i] >= keepAtOrAbove ? scores[i] : margins[i]).ThenByDescending(i => scores[i]).ThenBy(i => i).ToArray();
        double[] values = scores.OrderDescending().ToArray();
        var stacked = new double[n];
        for (int rank = 0; rank < n; rank++)
            stacked[order[rank]] = values[rank];
        return stacked;
    }

    /// <summary>
    /// Measures the rows near the threshold and re-orders them. Rows below the provisional floor keep their scores, and every value above it stays above it. So do rows at or under <see cref="KeptQValue"/>.
    /// </summary>
    /// <returns>The matches in their input order, with their scores reassigned.</returns>
    public static List<DiaPrecursorMatch> Apply(IReadOnlyList<DiaPrecursorMatch> matches, Func<int, PrecursorIons?> ionsOf, IReadOnlyList<(double Low, double High)> windowBounds,
        IReadOnlyList<MsDataScan[]> windowScans, double rtToleranceMinutes, double tolerancePpm, int apexHalfWidthScans, bool mostIntensePeak, int topFragments, int seed, int threads = 1, Action<string>? report = null)
    {
        double[]?[] columns = MeasureAll(matches, ionsOf, windowBounds, windowScans, rtToleranceMinutes, tolerancePpm, apexHalfWidthScans, mostIntensePeak, topFragments, threads);
        int[] measured = Enumerable.Range(0, matches.Count).Where(i => columns[i] is not null).ToArray();
        var applied = matches.ToList();
        if (measured.Length == 0)
            return applied;
        double[] q = TargetDecoyQValues.Compute(matches.Select(m => m.Score).ToArray(), matches.Select(m => m.IsDecoy).ToArray());
        double keepAtOrAbove = matches.Where((_, i) => q[i] <= KeptQValue).Select(m => m.Score).DefaultIfEmpty(double.PositiveInfinity).Min();
        double[] stacked = Stack(measured.Select(i => matches[i].Score).ToList(), measured.Select(i => columns[i]!).ToList(), measured.Select(i => matches[i].IsDecoy).ToList(),
            measured.Select(i => matches[i].FullSequence).ToList(), seed, threads, keepAtOrAbove);
        for (int k = 0; k < measured.Length; k++)
            applied[measured[k]] = matches[measured[k]] with { Score = stacked[k] };
        report?.Invoke($"neighbour evidence: {measured.Length} of {matches.Count} precursors measured and re-ranked ({measured.Count(i => matches[i].IsDecoy)} decoys)");
        return applied;
    }

    /// <summary>The square-root cosine of a spectrum against the library; 0 when either is empty.</summary>
    private static double Cosine(IReadOnlyList<(double Mz, double Intensity)> ions, double[] observed)
    {
        double dot = 0, normLibrary = 0, normObserved = 0;
        for (int i = 0; i < ions.Count; i++)
        {
            dot += Math.Sqrt(observed[i] * ions[i].Intensity);
            normLibrary += ions[i].Intensity;
            normObserved += observed[i];
        }
        return normLibrary > 0 && normObserved > 0 ? dot / Math.Sqrt(normLibrary * normObserved) : 0;
    }

    private static int FirstAtOrAfter(IReadOnlyList<Neighbour> neighbours, double rtMinutes)
    {
        int low = 0, high = neighbours.Count;
        while (low < high)
        {
            int middle = (low + high) / 2;
            if (neighbours[middle].ApexRt < rtMinutes) low = middle + 1; else high = middle;
        }
        return low;
    }

    private static int NearestScan(MsDataScan[] scans, double rtMinutes)
    {
        int low = 0, high = scans.Length;
        while (low < high)
        {
            int middle = (low + high) / 2;
            if (scans[middle].RetentionTime < rtMinutes) low = middle + 1; else high = middle;
        }
        if (low == 0)
            return 0;
        if (low == scans.Length)
            return scans.Length - 1;
        return rtMinutes - scans[low - 1].RetentionTime <= scans[low].RetentionTime - rtMinutes ? low - 1 : low;
    }
}
