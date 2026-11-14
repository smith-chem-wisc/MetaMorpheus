#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>
/// A run's MS1 mass offset (observed minus library m/z, in ppm) as a function of retention time: the median signed error in
/// equal-count RT bins, linearly interpolated between bin centres and held flat beyond the first and last. Fitted from
/// confident identifications (<see cref="DiaIrtSelfCalibration"/>) and applied by shifting every MS1 lookup
/// (<see cref="DiaLibrarySearchParameters.Ms1Offset"/>). On PXD022589 the offset drifts from +0.1 to +2.9 ppm over the
/// gradient, so a fixed window centred at 0 clips real precursors unevenly.
/// </summary>
public sealed class Ms1OffsetModel
{
    /// <summary>Errors per bin; fewer points give fewer bins, and under this many in total a single median.</summary>
    public const int PointsPerBin = 100;

    /// <summary>Most bins fitted.</summary>
    public const int MaximumBins = 20;

    private readonly double[] _rts;
    private readonly double[] _offsets;

    private Ms1OffsetModel(double[] rts, double[] offsets, double residualSpreadPpm)
    {
        _rts = rts;
        _offsets = offsets;
        ResidualSpreadPpm = residualSpreadPpm;
    }

    /// <summary>
    /// The errors' spread around the fitted offset, in ppm: 1.4826 x the median absolute residual (an SD for normal errors),
    /// what remains once the drift is removed. NaN for a constant model or one fitted from nothing.
    /// </summary>
    public double ResidualSpreadPpm { get; }

    /// <summary>The same offset at every retention time.</summary>
    public static Ms1OffsetModel Constant(double ppm) => new([0], [ppm], double.NaN);

    /// <summary>Fits the offset from (retention time in minutes, signed MS1 error in ppm) pairs; non-finite errors are ignored.</summary>
    public static Ms1OffsetModel Fit(IEnumerable<(double RtMinutes, double ErrorPpm)> points)
    {
        ArgumentNullException.ThrowIfNull(points);
        var sorted = points.Where(p => double.IsFinite(p.RtMinutes) && double.IsFinite(p.ErrorPpm)).OrderBy(p => p.RtMinutes).ToArray();
        if (sorted.Length == 0)
            return Constant(0);
        int bins = Math.Clamp(sorted.Length / PointsPerBin, 1, MaximumBins);
        var rts = new double[bins];
        var offsets = new double[bins];
        for (int b = 0; b < bins; b++)
        {
            int from = b * sorted.Length / bins, to = (b + 1) * sorted.Length / bins;
            var bin = sorted[from..to];
            rts[b] = Median(bin.Select(p => p.RtMinutes));
            offsets[b] = Median(bin.Select(p => p.ErrorPpm));
        }
        var fitted = new Ms1OffsetModel(rts, offsets, double.NaN);
        return new Ms1OffsetModel(rts, offsets, 1.4826 * Median(sorted.Select(p => Math.Abs(p.ErrorPpm - fitted.OffsetPpm(p.RtMinutes)))));
    }

    /// <summary>The offset at <paramref name="rtMinutes"/>, in ppm.</summary>
    public double OffsetPpm(double rtMinutes)
    {
        if (_rts.Length == 1 || rtMinutes <= _rts[0])
            return _offsets[0];
        if (rtMinutes >= _rts[^1])
            return _offsets[^1];
        int i = 1;
        while (_rts[i] < rtMinutes)
            i++;
        double span = _rts[i] - _rts[i - 1];
        return span <= 0 ? _offsets[i] : _offsets[i - 1] + (_offsets[i] - _offsets[i - 1]) * (rtMinutes - _rts[i - 1]) / span;
    }

    private static double Median(IEnumerable<double> values)
    {
        var sorted = values.Order().ToArray();
        int n = sorted.Length;
        return n % 2 == 1 ? sorted[n / 2] : (sorted[n / 2 - 1] + sorted[n / 2]) / 2;
    }
}
