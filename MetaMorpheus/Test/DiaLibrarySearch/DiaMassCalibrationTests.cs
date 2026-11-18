#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DiaLibrarySearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using Readers.SpectralLibrary;

namespace Test.DiaLibrarySearch;

/// <summary>
/// Mass calibration for DIA, by MetaMorpheus's own CalibrationEngine (the DDA calibration task's): the run's confident
/// identifications supply labelled points (fragments at the apex scan, the precursor's M0 in the nearest MS1 scan), and the
/// engine corrects every scan's m/z by its smoothed local error. Instruments other than PXD005573's drift by several ppm.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaMassCalibrationTests
{
    private string _directory = "";

    [OneTimeSetUp]
    public void SetUp()
    {
        _directory = Path.Combine(TestContext.CurrentContext.TestDirectory, "DiaMassCalibrationTests");
        Directory.CreateDirectory(_directory);
    }

    [OneTimeTearDown]
    public void TearDown() => Directory.Delete(_directory, true);

    private static double MedianPlantedFragmentPpm(IEnumerable<MsDataScan> scans, SyntheticDiaRun run)
    {
        var planted = run.Library.Where(e => run.PlantedSequences.Contains(e.FullSequence)).ToList();
        var errors = new List<double>();
        foreach (var scan in scans.Where(s => s.MsnOrder == 2))
            foreach (var entry in planted.Where(e => e.PrecursorMz >= scan.IsolationRange.Minimum && e.PrecursorMz < scan.IsolationRange.Maximum
                         && Math.Abs(scan.RetentionTime - SyntheticDiaRun.TrueRtMinutes(e.RetentionTime)) < 0.02))
                foreach (var f in entry.MatchedFragmentIons)
                {
                    int? i = scan.MassSpectrum.GetClosestPeakIndex(f.Mz);
                    if (i is int k)
                    {
                        double ppm = (scan.MassSpectrum.XArray[k] - f.Mz) / f.Mz * 1e6;
                        if (Math.Abs(ppm) < 20)
                            errors.Add(ppm);
                    }
                }
        return errors.Order().ElementAt(errors.Count / 2);
    }

    [Test]
    public void ARunOffsetBy8PpmIsCalibratedFromItsOwnIdentifications()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, withMs1: true, ppmOffset: 8);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var rt = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters());
        var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, rt.Model,
            new DiaLibrarySearchParameters(IrtHalfWindow: rt.IrtHalfWindow), new CommonParameters(), [], []).Run();
        Assert.That(MedianPlantedFragmentPpm(run.Scans, run), Is.EqualTo(8).Within(0.5), "the fixture's offset");

        var calibration = DiaMassCalibration.Calibrate(run.Scans, results.Matches, library, new CommonParameters());

        Assert.That(calibration.Ms2Points, Is.GreaterThan(100));
        Assert.That(calibration.Ms1Points, Is.GreaterThan(50));
        Assert.That(calibration.Ms2MedianPpmBefore, Is.EqualTo(8).Within(1));
        Assert.That(calibration.Ms1MedianPpmBefore, Is.EqualTo(8).Within(1));
        Assert.That(MedianPlantedFragmentPpm(calibration.Scans, run), Is.EqualTo(0).Within(1), "corrected");
        Assert.That(calibration.Scans, Has.Length.EqualTo(run.Scans.Length));
    }

    /// <summary>
    /// MetaMorpheus does not calibrate DIA, and its CalibrationEngine shows why: it also shifts each MS2 scan's isolation m/z
    /// by that scan's own correction, which for DDA is a measured precursor m/z. In DIA it is the instrument's fixed window,
    /// and shifting it per scan splits every window into one-scan windows, so nothing can be extracted. Only the spectra may
    /// change: the window, precursor fields and everything else stay as recorded, and the calibrated run still searches.
    /// </summary>
    [Test]
    public void CalibrationKeepsEveryIsolationWindowAndTheRunStillSearches()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, withMs1: true, ppmOffset: 8);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var rt = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters());
        var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, rt.Model,
            new DiaLibrarySearchParameters(IrtHalfWindow: rt.IrtHalfWindow), new CommonParameters(), [], []).Run();

        var calibrated = DiaMassCalibration.Calibrate(run.Scans, results.Matches, library, new CommonParameters()).Scans;

        for (int i = 0; i < run.Scans.Length; i++)
        {
            Assert.That(calibrated[i].IsolationMz, Is.EqualTo(run.Scans[i].IsolationMz), $"scan {i} isolation m/z");
            Assert.That(calibrated[i].IsolationWidth, Is.EqualTo(run.Scans[i].IsolationWidth));
            Assert.That(calibrated[i].SelectedIonMZ, Is.EqualTo(run.Scans[i].SelectedIonMZ));
            Assert.That(calibrated[i].RetentionTime, Is.EqualTo(run.Scans[i].RetentionTime));
            Assert.That(calibrated[i].OneBasedScanNumber, Is.EqualTo(run.Scans[i].OneBasedScanNumber));
        }
        var again = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(calibrated, library, rt.Model,
            new DiaLibrarySearchParameters(IrtHalfWindow: rt.IrtHalfWindow), new CommonParameters(), [], []).Run();
        int found = again.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence));
        Assert.That(found, Is.GreaterThanOrEqualTo((int)Math.Ceiling(0.9 * run.PlantedSequences.Count)));
    }

    private (string LibraryPath, MsDataScan[] Scans, List<DiaPrecursorMatch> Confident)? _searched;

    /// <summary>The 8 ppm run, searched once, and its targets at <see cref="DiaMassCalibration.LabelQValue"/>.</summary>
    private (string LibraryPath, MsDataScan[] Scans, List<DiaPrecursorMatch> Confident) Searched()
    {
        if (_searched is { } searched)
            return searched;
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, withMs1: true, ppmOffset: 8);
        string path = run.WriteLibrary(_directory);
        using var library = MslLibrary.Load(path);
        var rt = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters());
        var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, rt.Model,
            new DiaLibrarySearchParameters(IrtHalfWindow: rt.IrtHalfWindow), new CommonParameters(), [], []).Run();
        var confident = results.Matches.Where(m => !m.IsDecoy && m.QValue <= DiaMassCalibration.LabelQValue).ToList();
        _searched = (path, run.Scans, confident);
        return _searched.Value;
    }

    private static void AssertSamePoints(DiaMassCalibrationResult actual, DiaMassCalibrationResult expected)
    {
        Assert.That(expected.Ms2Points, Is.GreaterThan(100));
        Assert.That(actual.Ms1Points, Is.EqualTo(expected.Ms1Points));
        Assert.That(actual.Ms2Points, Is.EqualTo(expected.Ms2Points));
        Assert.That(actual.Ms1MedianPpmBefore, Is.EqualTo(expected.Ms1MedianPpmBefore));
        Assert.That(actual.Ms2MedianPpmBefore, Is.EqualTo(expected.Ms2MedianPpmBefore));
    }

    /// <summary>
    /// An identification labels from the scans nearest its apex time, which need not be a scan's own time: a time between two
    /// scans takes the nearer, and a time before the first scan or after the last takes that end's scan. Moving every apex a
    /// third of a cycle either way, and two of them beyond the run's ends, labels exactly the same points.
    /// </summary>
    [Test]
    public void AnApexBetweenOrBeyondTheScansLabelsFromTheNearestScan()
    {
        var (path, scans, confident) = Searched();
        using var library = MslLibrary.Load(path);
        double first = scans.Min(s => s.RetentionTime), last = scans.Max(s => s.RetentionTime);
        double shift = 0.3 * SyntheticDiaRun.CycleMinutes;
        var onScans = confident
            .Append(confident[0] with { ApexRt = new RtMinutes(first) })
            .Append(confident[0] with { ApexRt = new RtMinutes(last) }).ToList();
        var offScans = confident.Select((m, i) => m with { ApexRt = new RtMinutes(m.ApexRt.Value + (i % 2 == 0 ? shift : -shift)) })
            .Append(confident[0] with { ApexRt = new RtMinutes(first - 1) })
            .Append(confident[0] with { ApexRt = new RtMinutes(last + 1) }).ToList();

        var expected = DiaMassCalibration.Calibrate(scans, onScans, library, new CommonParameters());
        var actual = DiaMassCalibration.Calibrate(scans, offScans, library, new CommonParameters());

        AssertSamePoints(actual, expected);
    }

    /// <summary>
    /// Only targets at <see cref="DiaMassCalibration.LabelQValue"/> label points. Decoys, targets above it, and a precursor that
    /// no isolation window holds and no MS1 peak lies near add nothing.
    /// </summary>
    [Test]
    public void OnlyConfidentTargetsLabelPoints()
    {
        var (path, scans, confident) = Searched();
        using var library = MslLibrary.Load(path);
        var unusable = confident.Select(m => m with { IsDecoy = true })
            .Concat(confident.Select(m => m with { QValue = 0.5 }))
            .Append(confident[0] with { PrecursorMz = 100 });

        var expected = DiaMassCalibration.Calibrate(scans, confident, library, new CommonParameters());
        var actual = DiaMassCalibration.Calibrate(scans, confident.Concat(unusable).ToList(), library, new CommonParameters());

        AssertSamePoints(actual, expected);
    }

    /// <summary>With no confident identifications there is nothing to calibrate on: the scans come back unchanged.</summary>
    [Test]
    public void WithoutIdentificationsTheScansAreReturnedUnchanged()
    {
        var run = SyntheticDiaRun.Build(50, _ => false, withMs1: true, ppmOffset: 8);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaMassCalibration.Calibrate(run.Scans, [], library, new CommonParameters());

        Assert.That(calibration.Scans, Is.SameAs(run.Scans));
        Assert.That(calibration.Ms2Points, Is.EqualTo(0));
    }
}
