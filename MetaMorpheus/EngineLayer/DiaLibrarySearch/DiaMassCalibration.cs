#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using EngineLayer.Calibration;
using MassSpectrometry;
using Readers;
using Readers.SpectralLibrary;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>A DIA run's scans after mass calibration, and the errors it was fitted on.</summary>
/// <param name="Scans">The calibrated scans, or the input array itself when there was nothing to calibrate on.</param>
/// <param name="Ms1MedianPpmBefore">Median precursor (M0) error of the labelled points, before correction.</param>
/// <param name="Ms2MedianPpmBefore">Median fragment error of the labelled points, before correction.</param>
public sealed record DiaMassCalibrationResult(MsDataScan[] Scans, int Ms1Points, int Ms2Points, double Ms1MedianPpmBefore, double Ms2MedianPpmBefore);

/// <summary>
/// Mass calibration for a DIA run by MetaMorpheus's DDA <see cref="CalibrationEngine"/>, unchanged. This class only turns DIA
/// identifications into its labelled points: for each target at <see cref="LabelQValue"/>, the six most intense library fragments
/// found within <see cref="SearchPpm"/> in the apex MS2 scan, and the precursor's monoisotopic peak in the nearest MS1 scan. The
/// engine averages each scan's errors (intensity-weighted), fills scans without points from their neighbours, smooths over
/// neighbouring scans, and shifts every m/z by the result.
/// </summary>
public static class DiaMassCalibration
{
    /// <summary>Identifications at or below this q-value label the calibration.</summary>
    public const double LabelQValue = 0.01;

    /// <summary>How far from the library m/z a peak may be and still label a point: wide enough for an uncalibrated run.</summary>
    public const double SearchPpm = 20;

    private const int FragmentsPerIdentification = 6;

    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public static DiaMassCalibrationResult Calibrate(MsDataScan[] scans, IReadOnlyList<DiaPrecursorMatch> matches, MslLibrary library,
        CommonParameters commonParameters)
    {
        ArgumentNullException.ThrowIfNull(scans);
        ArgumentNullException.ThrowIfNull(matches);
        ArgumentNullException.ThrowIfNull(library);
        ArgumentNullException.ThrowIfNull(commonParameters);

        var ms1 = scans.Where(s => s.MsnOrder == 1).OrderBy(s => s.RetentionTime).ToArray();
        double[] ms1Rts = ms1.Select(s => s.RetentionTime).ToArray();
        var ms2ByWindow = scans.Where(s => s.MsnOrder == 2 && s.IsolationRange is not null)
            .GroupBy(s => (s.IsolationRange.Minimum, s.IsolationRange.Maximum))
            .Select(g => (Low: g.Key.Minimum, High: g.Key.Maximum, Scans: g.OrderBy(s => s.RetentionTime).ToArray()))
            .ToArray();

        var ms1Points = new List<LabeledDataPoint>();
        var ms2Points = new List<LabeledDataPoint>();
        foreach (var match in matches.Where(m => !m.IsDecoy && m.QValue <= LabelQValue))
        {
            var entry = library.GetEntry(match.PrecursorIndex);
            if (entry is null)
                continue;
            var window = ms2ByWindow.FirstOrDefault(w => match.PrecursorMz >= w.Low && match.PrecursorMz < w.High);
            if (window.Scans is { Length: > 0 })
            {
                var apex = Nearest(window.Scans, s => s.RetentionTime, match.ApexRt.Value);
                foreach (var fragment in entry.MatchedFragmentIons.OrderByDescending(f => f.Intensity).Take(FragmentsPerIdentification))
                    if (Label(apex, fragment.Mz) is { } point)
                        ms2Points.Add(point);
            }
            if (ms1.Length > 0)
            {
                var survey = ms1[NearestIndex(ms1Rts, match.ApexRt.Value)];
                if (Label(survey, match.PrecursorMz) is { } point)
                    ms1Points.Add(point);
            }
        }
        if (ms1Points.Count == 0 || ms2Points.Count == 0)
            return new DiaMassCalibrationResult(scans, ms1Points.Count, ms2Points.Count, double.NaN, double.NaN);

        // The engine reads each list in scan order
        ms1Points.Sort((a, b) => a.ScanNumber.CompareTo(b.ScanNumber));
        ms2Points.Sort((a, b) => a.ScanNumber.CompareTo(b.ScanNumber));
        var datapoints = new DataPointAquisitionResults(null!, [], ms1Points, ms2Points, 0, 0, 0, 0, commonParameters);
        var file = new GenericMsDataFile(scans, new SourceFile("no nativeID format", "mzML format", null, null, null));
        var engine = new CalibrationEngine(file, datapoints, commonParameters, [], []);
        engine.Run();
        // Keep only the engine's corrected spectra. It also shifts each MS2 scan's isolation and selected-ion m/z, which for
        // DDA are measured precursor m/z; in DIA they are the instrument's fixed window, and shifting them per scan would
        // split every window into one-scan windows (MetaMorpheus does not calibrate DIA)
        var corrected = engine.CalibratedDataFile.GetAllScansList();
        var calibrated = new MsDataScan[scans.Length];
        for (int i = 0; i < scans.Length; i++)
            calibrated[i] = WithSpectrum(scans[i], corrected[i].MassSpectrum);
        return new DiaMassCalibrationResult(calibrated, ms1Points.Count, ms2Points.Count,
            Median(ms1Points.Select(p => p.RelativeMzError * 1e6)), Median(ms2Points.Select(p => p.RelativeMzError * 1e6)));
    }

    /// <summary>A labelled point from the peak nearest <paramref name="mz"/>, if it lies within <see cref="SearchPpm"/>; otherwise null.</summary>
    private static LabeledDataPoint? Label(MsDataScan scan, double mz)
    {
        int? index = scan.MassSpectrum.GetClosestPeakIndex(mz);
        if (index is not int i || Math.Abs(scan.MassSpectrum.XArray[i] - mz) > mz * SearchPpm * 1e-6)
            return null;
        double intensity = scan.MassSpectrum.YArray[i];
        return new LabeledDataPoint(scan.MassSpectrum.XArray[i], scan.OneBasedScanNumber, Math.Log10(Math.Max(1, scan.TotalIonCurrent)),
            Math.Log10(Math.Max(1e-3, scan.InjectionTime ?? 1)), Math.Log10(Math.Max(1, intensity)), mz, null);
    }

    /// <summary>The scan as recorded, with a different spectrum.</summary>
    private static MsDataScan WithSpectrum(MsDataScan scan, MzSpectrum spectrum) => new(
        spectrum, scan.OneBasedScanNumber, scan.MsnOrder, scan.IsCentroid, scan.Polarity, scan.RetentionTime, scan.ScanWindowRange,
        scan.ScanFilter, scan.MzAnalyzer, scan.TotalIonCurrent, scan.InjectionTime, scan.NoiseData, scan.NativeId, scan.SelectedIonMZ,
        scan.SelectedIonChargeStateGuess, scan.SelectedIonIntensity, scan.IsolationMz, scan.IsolationWidth, scan.DissociationType,
        scan.OneBasedPrecursorScanNumber, scan.SelectedIonMonoisotopicGuessMz);

    private static MsDataScan Nearest(MsDataScan[] sorted, Func<MsDataScan, double> key, double value) =>
        sorted[NearestIndex(sorted.Select(key).ToArray(), value)];

    private static int NearestIndex(double[] sorted, double value)
    {
        int i = Array.BinarySearch(sorted, value);
        if (i >= 0)
            return i;
        i = ~i;
        if (i == 0)
            return 0;
        if (i == sorted.Length)
            return sorted.Length - 1;
        return value - sorted[i - 1] <= sorted[i] - value ? i - 1 : i;
    }

    private static double Median(IEnumerable<double> values)
    {
        var sorted = values.Order().ToArray();
        return sorted.Length == 0 ? double.NaN : sorted[sorted.Length / 2];
    }
}
