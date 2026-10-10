#nullable enable
using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;

namespace Test.DiaLibrarySearch;

/// <summary>
/// The synthetic run has to be right before it can judge the engine: planted signal where the truth says, nowhere
/// else, and in the window that isolates it.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class SyntheticDiaRunTests
{
    [Test]
    public void EveryTargetHasAReversedDecoyAtTheSamePrecursorAndIrt()
    {
        var run = SyntheticDiaRun.Build(20, _ => false);

        var targets = run.Library.Where(e => !e.IsDecoy).ToList();
        var decoys = run.Library.Where(e => e.IsDecoy).ToList();
        Assert.That(targets.Count, Is.EqualTo(20));
        Assert.That(decoys.Count, Is.EqualTo(20));
        foreach (var target in targets)
        {
            var decoy = decoys.Single(d => d.FullSequence == new string(target.FullSequence.Reverse().ToArray()));
            Assert.That(decoy.FullSequence, Is.Not.EqualTo(target.FullSequence));
            Assert.That(decoy.PrecursorMz, Is.EqualTo(target.PrecursorMz));
            Assert.That(decoy.RetentionTime, Is.EqualTo(target.RetentionTime));
        }
    }

    [Test]
    public void EveryScanIsolatesOneWindowOfTheScheme()
    {
        var run = SyntheticDiaRun.Build(5, _ => false);

        var windows = run.Scans.Select(s => (s.IsolationRange.Minimum, s.IsolationRange.Maximum)).Distinct().ToList();
        Assert.That(windows.Count, Is.EqualTo(SyntheticDiaRun.WindowCount));
        Assert.That(run.Scans.All(s => s.MsnOrder == 2));
        Assert.That(windows.Min(w => w.Minimum), Is.EqualTo(SyntheticDiaRun.FirstWindowLowMz).Within(1e-9));
    }

    /// <summary>At the true apex, the scan whose window holds the precursor carries each fragment at its library ratio.</summary>
    [Test]
    public void APlantedPrecursorElutesAtItsTrueRetentionTimeInItsOwnWindow()
    {
        var run = SyntheticDiaRun.Build(30, e => !e.IsDecoy, noisePeaksPerScan: 0);
        var entry = run.Library.First(e => !e.IsDecoy);
        double apexRt = SyntheticDiaRun.TrueRtMinutes(entry.RetentionTime);
        var tolerance = new PpmTolerance(5);

        var apexScan = run.Scans
            .Where(s => s.IsolationRange.Contains(entry.PrecursorMz))
            .OrderBy(s => Math.Abs(s.RetentionTime - apexRt))
            .First();
        double Intensity(float mz)
        {
            int i = apexScan.MassSpectrum.GetClosestPeakIndex(mz);
            return tolerance.Within(apexScan.MassSpectrum.XArray[i], mz) ? apexScan.MassSpectrum.YArray[i] : 0;
        }

        var observed = entry.MatchedFragmentIons.Select(f => Intensity(f.Mz)).ToArray();
        var expected = entry.MatchedFragmentIons.Select(f => (double)f.Intensity).ToArray();
        Assert.That(observed.All(i => i > 0));
        for (int f = 1; f < observed.Length; f++)
            Assert.That(observed[f] / observed[0], Is.EqualTo(expected[f] / expected[0]).Within(0.02), $"fragment {f} ratio");

        var otherWindow = run.Scans.First(s => !s.IsolationRange.Contains(entry.PrecursorMz)
            && Math.Abs(s.RetentionTime - apexScan.RetentionTime) < 1e-9);
        // With noise off, a window isolating nothing planted is empty (GetClosestPeakIndex does not guard that)
        var other = otherWindow.MassSpectrum;
        bool fragmentThere = other.Size > 0
            && tolerance.Within(other.XArray[other.GetClosestPeakIndex(entry.MatchedFragmentIons[0].Mz)], entry.MatchedFragmentIons[0].Mz);
        Assert.That(fragmentThere, Is.False, "the fragment must not appear in a window that does not isolate its precursor");
    }

    [Test]
    public void TheTruthMapIsNonlinearAndIncreasingOverTheLibraryRange()
    {
        double previous = double.NegativeInfinity;
        for (double irt = -20; irt <= 120; irt += 1)
        {
            double rt = SyntheticDiaRun.TrueRtMinutes(irt);
            Assert.That(rt, Is.GreaterThan(previous));
            Assert.That(rt, Is.InRange(0, SyntheticDiaRun.RunMinutes));
            previous = rt;
        }
        double slopeLow = SyntheticDiaRun.TrueRtMinutes(-19) - SyntheticDiaRun.TrueRtMinutes(-20);
        double slopeHigh = SyntheticDiaRun.TrueRtMinutes(120) - SyntheticDiaRun.TrueRtMinutes(119);
        Assert.That(slopeHigh / slopeLow, Is.GreaterThan(1.5), "a linear truth would let iRT-as-minutes bugs hide");
    }

    [Test]
    public void TheSameSeedBuildsTheSameRun()
    {
        var a = SyntheticDiaRun.Build(10, e => !e.IsDecoy, seed: 7);
        var b = SyntheticDiaRun.Build(10, e => !e.IsDecoy, seed: 7);

        Assert.That(a.Library.Select(e => e.FullSequence), Is.EqualTo(b.Library.Select(e => e.FullSequence)));
        Assert.That(a.Scans[123].MassSpectrum.YArray, Is.EqualTo(b.Scans[123].MassSpectrum.YArray));
    }
}
