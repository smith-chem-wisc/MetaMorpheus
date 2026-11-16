#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>
/// What a run's confident targets say of the library's predicted fragment intensities, by ion type, distance from the terminus and ion charge. A predicted library is fixed at one collision energy on
/// some other instrument, and on both benchmark runs its error has a shape: y3 reads about twice the prediction, falling steadily to about 0.6x by y14. A table fitted on half of the confident
/// targets raises the other half's 10th-percentile cosine by 0.03 (PXD022589 HF-X) and 0.06 (PXD005573 1 h). Calibration fits it, and the search then reads the library through it, in choosing a
/// precursor's strongest fragments as well as in scoring them.
/// </summary>
/// <remarks>
/// Only targets fit the table, but it goes by ion class alone, so a decoy's ions are moved exactly as a target's are.
/// </remarks>
public sealed class DiaFragmentIntensityModel
{
    /// <summary>The furthest a multiplier may move an ion, either way.</summary>
    public const double MaximumMultiplier = 4;

    /// <summary>A weaker library ion does not fit the table, since it is only counted when it happens to be seen, which biases its class upward.</summary>
    public const double LibraryFloor = 0.05;

    /// <summary>The fewest ions a class needs before it is moved.</summary>
    public const int MinimumIonsPerClass = 30;

    /// <summary>Ions at least this far from their terminus share a class.</summary>
    private const int FurthestNumber = 14;

    /// <summary>The fragment rule's bounds (<see cref="DiaLibrarySearchParameters.FragmentMinimumResidues"/>, <see cref="DiaLibrarySearchParameters.FragmentMinimumMz"/>): nearer ions are never scored by default.</summary>
    private const int NearestNumber = 3;

    private const double LowestMz = 200;

    private readonly Dictionary<(ProductType Type, int Number, int Charge), double> _multipliers;

    private DiaFragmentIntensityModel(Dictionary<(ProductType, int, int), double> multipliers) => _multipliers = multipliers;

    /// <summary>What the table does to an ion's library intensity; 1 for an ion outside the table or a class seen too rarely.</summary>
    public double Multiplier(MslFragmentIon ion) =>
        InTheTable(ion) && _multipliers.TryGetValue(ClassOf(ion), out double multiplier) ? multiplier : 1;

    /// <summary>An ion's library intensity, read through the table.</summary>
    public double Intensity(MslFragmentIon ion) => ion.Intensity * Multiplier(ion);

    /// <summary>Fits the table on confident targets' ions: per class, the median ratio of an ion's share of what its precursor shows to its share of what the library predicts, both among the ions that fit.</summary>
    /// <param name="anchors">Each confident target's library ions, with what its apex holds at each.</param>
    public static DiaFragmentIntensityModel Fit(IEnumerable<IReadOnlyList<(MslFragmentIon Ion, double Observed)>> anchors)
    {
        ArgumentNullException.ThrowIfNull(anchors);
        var logRatios = new Dictionary<(ProductType, int, int), List<double>>();
        foreach (var anchor in anchors)
        {
            var seen = anchor.Where(read => InTheTable(read.Ion) && read.Ion.Intensity >= LibraryFloor && read.Observed > 0).ToList();
            double observedTotal = seen.Sum(read => read.Observed);
            double libraryTotal = seen.Sum(read => (double)read.Ion.Intensity);
            foreach (var (ion, observed) in seen)
            {
                var key = ClassOf(ion);
                if (!logRatios.TryGetValue(key, out var ratios))
                    logRatios[key] = ratios = [];
                ratios.Add(Math.Log(observed / observedTotal / (ion.Intensity / libraryTotal)));
            }
        }
        double limit = Math.Log(MaximumMultiplier);
        return new DiaFragmentIntensityModel(logRatios.Where(c => c.Value.Count >= MinimumIonsPerClass)
            .ToDictionary(c => c.Key, c => Math.Exp(Math.Clamp(c.Value.Order().ElementAt(c.Value.Count / 2), -limit, limit))));
    }

    /// <summary>What a confident target's apex holds at each of its library ions, over three scans, in the window whose scans come closest to the apex.</summary>
    /// <param name="entry">The target.</param>
    /// <param name="apexRt">Its apex.</param>
    /// <param name="windowsHoldingIt">The scans of each isolation window holding its m/z, each in retention-time order.</param>
    /// <param name="tolerancePpm">The fragment tolerance.</param>
    /// <param name="mostIntensePeak">Take the most intense peak in tolerance rather than the closest, as the search does.</param>
    /// <returns>Its ions with what was read; empty with no window to read.</returns>
    public static List<(MslFragmentIon Ion, double Observed)> ReadAnchor(MslLibraryEntry entry, RtMinutes apexRt, IReadOnlyList<MsDataScan[]> windowsHoldingIt, double tolerancePpm, bool mostIntensePeak)
    {
        MsDataScan[]? window = null;
        int apex = -1;
        double gap = double.PositiveInfinity;
        foreach (var scans in windowsHoldingIt)
        {
            if (scans.Length == 0)
                continue;
            int scan = NearestScan(scans, apexRt.Value);
            double distance = Math.Abs(scans[scan].RetentionTime - apexRt.Value);
            if (distance < gap)
            {
                gap = distance;
                window = scans;
                apex = scan;
            }
        }
        if (window is null)
            return [];
        var tolerance = new PpmTolerance(tolerancePpm);
        var read = new List<(MslFragmentIon, double)>(entry.MatchedFragmentIons.Count);
        foreach (var ion in entry.MatchedFragmentIons)
        {
            double sum = 0;
            for (int k = Math.Max(0, apex - 1); k <= Math.Min(window.Length - 1, apex + 1); k++)
            {
                var spectrum = window[k].MassSpectrum;
                if (spectrum.Size == 0)
                    continue;
                int peak = DiaLibrarySearchEngine.FragmentPeakIndex(spectrum, ion.Mz, tolerance, mostIntensePeak);
                if (peak >= 0)
                    sum += spectrum.YArray[peak];
            }
            read.Add((ion, sum));
        }
        return read;
    }

    private static bool InTheTable(MslFragmentIon ion) =>
        !ion.IsInternalFragment && ion.NeutralLoss == 0 && ion.ProductType is ProductType.y or ProductType.b && ion.FragmentNumber >= NearestNumber && ion.Mz >= LowestMz;

    private static (ProductType, int, int) ClassOf(MslFragmentIon ion) => (ion.ProductType, Math.Min(ion.FragmentNumber, FurthestNumber), ion.Charge);

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
