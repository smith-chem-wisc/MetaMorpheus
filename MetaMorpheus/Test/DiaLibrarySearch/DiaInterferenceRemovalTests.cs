#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using EngineLayer.DiaLibrarySearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;

namespace Test.DiaLibrarySearch;

/// <summary>
/// Interference removal, after DIA-NN's remove_ifs: a match whose top library fragments a better-scoring match co-eluting
/// in the same or a neighbouring window already explains is dropped, target or decoy alike. On the whole-proteome HeLa
/// library such "riders" were most of the decoys scoring above the 2% target threshold.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaInterferenceRemovalTests
{
    private const double Ppm = 20;
    private const double RtTolerance = 0.1;
    private const int Explained = 4;

    // Precursor 0's fragments; the others copy some of them
    private static readonly float[] Base = [300f, 400f, 500f, 600f, 700f, 800f];

    private static float[] Sharing(int shared) => Base.Take(shared).Concat(Enumerable.Range(0, 6 - shared).Select(i => 1100f + 37 * i)).ToArray();

    private static DiaPrecursorMatch Match(int index, double score, bool isDecoy = false, double rt = 10, double mz = 500) =>
        new(index, $"PEPTIDE{index}K", 2, mz, isDecoy, new Irt(0), new RtMinutes(rt), new Irt(0), score, [score]);

    // 200 unrelated, well-separated targets keep the provisional floor low, so the cases under test take part
    private static List<DiaPrecursorMatch> Background() =>
        Enumerable.Range(100, 200).Select(i => Match(i, 1 + i / 1000.0, rt: 20 + i)).ToList();

    private static List<DiaPrecursorMatch> Remove(List<DiaPrecursorMatch> matches, Dictionary<int, float[]> fragments) =>
        DiaInterferenceRemoval.Remove(matches,
            index => fragments.TryGetValue(index, out var mzs) ? mzs : Enumerable.Range(0, 6).Select(i => 2000f + index * 7 + i).ToArray(),
            mz => (int)((mz - 400) / 25), RtTolerance, Ppm, Explained);

    [Test]
    public void ADecoyRidingOnABetterTargetIsRemovedAndTheTargetKept()
    {
        var matches = Background().Append(Match(0, 5)).Append(Match(1, 4.9, isDecoy: true)).ToList();

        var kept = Remove(matches, new() { [0] = Base, [1] = Sharing(5) });

        Assert.That(kept.Select(m => m.PrecursorIndex), Does.Contain(0).And.Not.Contain(1));
        Assert.That(kept, Has.Count.EqualTo(matches.Count - 1));
    }

    /// <summary>Removal is blind to the label: a target explained by a better decoy goes too.</summary>
    [Test]
    public void ATargetExplainedByABetterDecoyIsRemovedToo()
    {
        var matches = Background().Append(Match(0, 5, isDecoy: true)).Append(Match(1, 4.9)).ToList();

        var kept = Remove(matches, new() { [0] = Base, [1] = Sharing(4) });

        Assert.That(kept.Select(m => m.PrecursorIndex), Does.Contain(0).And.Not.Contain(1));
    }

    [TestCase(3, true, TestName = "Three shared fragments are coincidence")]
    [TestCase(4, false, TestName = "Four shared fragments are interference")]
    public void TheThresholdIsOnSharedTopFragments(int shared, bool survives)
    {
        var matches = Background().Append(Match(0, 5)).Append(Match(1, 4.9)).ToList();

        var kept = Remove(matches, new() { [0] = Base, [1] = Sharing(shared) });

        Assert.That(kept.Any(m => m.PrecursorIndex == 1), Is.EqualTo(survives));
    }

    [TestCase(10.05, 500, false, TestName = "Co-eluting in the same window")]
    [TestCase(10.5, 500, true, TestName = "Eluting apart")]
    [TestCase(10.05, 530, false, TestName = "In the neighbouring window")]
    [TestCase(10.05, 560, true, TestName = "Two windows away")]
    public void OnlyACoElutingMatchInTheSameOrANeighbouringWindowExplains(double rt, double mz, bool survives)
    {
        var matches = Background().Append(Match(0, 5)).Append(Match(1, 4.9, rt: rt, mz: mz)).ToList();

        var kept = Remove(matches, new() { [0] = Base, [1] = Sharing(6) });

        Assert.That(kept.Any(m => m.PrecursorIndex == 1), Is.EqualTo(survives));
    }

    /// <summary>The better match explains; the worse one, even when it copies the better, never removes it.</summary>
    [Test]
    public void TheBetterScoringMatchIsTheOneKept()
    {
        var matches = Background().Append(Match(0, 4.9)).Append(Match(1, 5)).ToList();

        var kept = Remove(matches, new() { [0] = Base, [1] = Sharing(6) });

        Assert.That(kept.Select(m => m.PrecursorIndex), Does.Contain(1).And.Not.Contain(0));
    }

    /// <summary>Matches too weak to be reported take no part, either as explainers or as explained.</summary>
    [Test]
    public void MatchesBelowTheProvisionalFloorAreLeftAlone()
    {
        var decoys = Enumerable.Range(400, 200).Select(i => Match(i, -5 + i / 1000.0, isDecoy: true, rt: 50 + i)).ToList();
        var matches = Background().Concat(decoys).Append(Match(0, -9)).Append(Match(1, -9.1)).ToList();

        var kept = Remove(matches, new() { [0] = Base, [1] = Sharing(6) });

        Assert.That(kept, Has.Count.EqualTo(matches.Count));
    }

    private static List<DiaPrecursorMatch> RemoveSameMzOnly(List<DiaPrecursorMatch> matches, Dictionary<int, float[]> fragments) =>
        DiaInterferenceRemoval.Remove(matches,
            index => fragments.TryGetValue(index, out var mzs) ? mzs : Enumerable.Range(0, 6).Select(i => 2000f + index * 7 + i).ToArray(),
            mz => (int)((mz - 400) / 25), RtTolerance, Ppm, Explained, sameMzOnly: true);

    /// <summary>
    /// DIA-NN removes a precursor only when the better one has its m/z, or has it as the +1 isotope (the lower one is then
    /// the better one's isotope peak). On the whole-proteome library our looser rule removed 239 of DIA-NN's identifications,
    /// 42% of which AlphaDIA also finds. With sameMzOnly, a different co-eluting precursor never removes another.
    /// </summary>
    [TestCase(500.0, false, TestName = "Same precursor m/z is removed")]
    [TestCase(500.0 + 1.0033548 / 2, false, TestName = "On the better precursor's +1 isotope is removed")]
    [TestCase(510.0, true, TestName = "A different m/z in the same window is kept")]
    public void WithSameMzOnlyADifferentPrecursorNeverRemovesAnother(double mz, bool survives)
    {
        var matches = Background().Append(Match(0, 5, mz: 500)).Append(Match(1, 4.9, mz: mz)).ToList();

        var kept = RemoveSameMzOnly(matches, new() { [0] = Base, [1] = Sharing(6) });

        Assert.That(kept.Any(m => m.PrecursorIndex == 1), Is.EqualTo(survives));
    }

    private static List<DiaPrecursorMatch> RemoveWithOwnIons(List<DiaPrecursorMatch> matches, Dictionary<int, float[]> fragments,
        Func<DiaPrecursorMatch, DiaPrecursorMatch, double> ownIons, double floor = 0.1, Action<string>? report = null) =>
        DiaInterferenceRemoval.Remove(matches,
            index => fragments.TryGetValue(index, out var mzs) ? mzs : Enumerable.Range(0, 6).Select(i => 2000f + index * 7 + i).ToArray(),
            mz => (int)((mz - 400) / 25), RtTolerance, Ppm, Explained, report, ownIonCoElution: ownIons, ownIonFloor: floor);

    /// <summary>
    /// Two library precursors of one m/z that share most fragments read the same signal, and the better score only says which fits it better. If the
    /// winner's own ions do not co-elute with what they share, nothing shows it is the real one. On the paired library two-thirds of the entrapments at
    /// 1% q are such winners, against 4% of targets. The rival stays removed: telling the two apart by these ions was worse than the score.
    /// </summary>
    [Test]
    public void AWinnerWhoseOwnIonsDoNotCoEluteGoesWithItsRival()
    {
        var matches = Background().Append(Match(0, 5)).Append(Match(1, 4.9)).ToList();

        var kept = RemoveWithOwnIons(matches, new() { [0] = Base, [1] = Sharing(5) }, (_, _) => 0.05);

        Assert.That(kept.Select(m => m.PrecursorIndex), Does.Not.Contain(0).And.Not.Contain(1));
        Assert.That(kept, Has.Count.EqualTo(matches.Count - 2));
    }

    [Test]
    public void AWinnerWhoseOwnIonsCoEluteStaysAndItsRivalGoes()
    {
        var matches = Background().Append(Match(0, 5)).Append(Match(1, 4.9)).ToList();

        var kept = RemoveWithOwnIons(matches, new() { [0] = Base, [1] = Sharing(5) }, (_, _) => 0.8);

        Assert.That(kept.Select(m => m.PrecursorIndex), Does.Contain(0).And.Not.Contain(1));
    }

    /// <summary>The rule asks the winner about its own ions, never the removed side about its.</summary>
    [Test]
    public void TheWinnerIsAskedAboutItsRival()
    {
        var matches = Background().Append(Match(0, 5)).Append(Match(1, 4.9)).ToList();
        var asked = new List<(int Winner, int Rival)>();

        RemoveWithOwnIons(matches, new() { [0] = Base, [1] = Sharing(5) }, (winner, rival) =>
        {
            asked.Add((winner.PrecursorIndex, rival.PrecursorIndex));
            return 0.8;
        });

        Assert.That(asked, Is.EqualTo(new[] { (0, 1) }));
    }

    /// <summary>Most targets have no rival. They have nothing to be told apart from, so they are never asked and their ions cost nothing.</summary>
    [Test]
    public void AMatchWithNoRivalIsNeverAsked()
    {
        var matches = Background().Append(Match(0, 5)).ToList();

        var kept = RemoveWithOwnIons(matches, new() { [0] = Base }, (_, _) => throw new InvalidOperationException("asked"));

        Assert.That(kept, Has.Count.EqualTo(matches.Count));
    }

    /// <summary>A winner must be told apart from every rival it beat.</summary>
    [Test]
    public void TheWeakestRivalDecides()
    {
        var matches = Background().Append(Match(0, 5)).Append(Match(1, 4.9)).Append(Match(2, 4.8)).ToList();

        var kept = RemoveWithOwnIons(matches, new() { [0] = Base, [1] = Sharing(5), [2] = Sharing(5) }, (_, rival) => rival.PrecursorIndex == 1 ? 0.8 : 0.0);

        Assert.That(kept.Select(m => m.PrecursorIndex), Does.Not.Contain(0));
    }

    /// <summary>Blind to the label, as removal is: a decoy winner goes the same way, or the target-decoy estimate would not stay fair.</summary>
    [Test]
    public void ADecoyWinnerGoesTheSameWay()
    {
        var matches = Background().Append(Match(0, 5, isDecoy: true)).Append(Match(1, 4.9)).ToList();

        var kept = RemoveWithOwnIons(matches, new() { [0] = Base, [1] = Sharing(5) }, (_, _) => 0.05);

        Assert.That(kept.Select(m => m.PrecursorIndex), Does.Not.Contain(0).And.Not.Contain(1));
    }

    /// <summary>When two better matches could explain it, the rival belongs to the one explaining the most fragments.</summary>
    [Test]
    public void TheRivalBelongsToTheWinnerExplainingTheMost()
    {
        float[] removedIons = Base;
        float[] explainsFour = Sharing(4);
        // Its top six are its own, so the first winner does not explain it, and the rest hold everything the removed match has
        float[] explainsAll = [.. Enumerable.Range(0, 6).Select(i => 3000f + 11 * i), .. Base];
        var matches = Background().Append(Match(0, 5)).Append(Match(1, 4.9)).Append(Match(2, 4.8)).ToList();
        var asked = new List<int>();

        RemoveWithOwnIons(matches, new() { [0] = explainsFour, [1] = explainsAll, [2] = removedIons }, (winner, rival) =>
        {
            if (rival.PrecursorIndex == 2)
                asked.Add(winner.PrecursorIndex);
            return 0.8;
        });

        Assert.That(asked, Is.EqualTo(new[] { 1 }));
    }

    /// <summary>A floor of 0 is the plain rule, and it never reads the ions.</summary>
    [Test]
    public void WithNoFloorNobodyIsAskedAndTheWinnerStays()
    {
        var matches = Background().Append(Match(0, 5)).Append(Match(1, 4.9)).ToList();

        var kept = RemoveWithOwnIons(matches, new() { [0] = Base, [1] = Sharing(5) }, (_, _) => throw new InvalidOperationException("asked"), floor: 0);

        Assert.That(kept.Select(m => m.PrecursorIndex), Does.Contain(0).And.Not.Contain(1));
    }

    /// <summary>Not a number means there was nothing to measure, which is no evidence against the winner.</summary>
    [Test]
    public void NothingToMeasureKeepsTheWinner()
    {
        var matches = Background().Append(Match(0, 5)).Append(Match(1, 4.9)).ToList();

        var kept = RemoveWithOwnIons(matches, new() { [0] = Base, [1] = Sharing(5) }, (_, _) => double.NaN);

        Assert.That(kept.Select(m => m.PrecursorIndex), Does.Contain(0).And.Not.Contain(1));
    }

    [Test]
    public void TheWinnersDroppedAreReported()
    {
        var matches = Background().Append(Match(0, 5, isDecoy: true)).Append(Match(1, 4.9)).Append(Match(2, 4.7, rt: 15)).Append(Match(3, 4.6, rt: 15)).ToList();
        var lines = new List<string>();

        RemoveWithOwnIons(matches, new() { [0] = Base, [1] = Sharing(5), [2] = Base, [3] = Sharing(5) }, (_, _) => 0.05, report: lines.Add);

        Assert.That(lines, Has.Some.Contains("2 winners dropped, own ions below 0.1 (1 decoys)"));
    }

    // A chromatographic peak of 21 scans with its apex at scan 10; the own ions below are placed on it, or away from it
    private const double PeakApexRt = 10;

    private static double Peak(int scan, int apex = 10) => 100_000 * Math.Exp(-(scan - apex) * (scan - apex) / 8.0);

    private static MsDataScan[] Window(params (double Mz, Func<int, double> Height)[] ions) =>
        Enumerable.Range(0, 21).Select(scan =>
        {
            var peaks = ions.Select(i => (i.Mz, Height: i.Height(scan))).Where(p => p.Height > 0).OrderBy(p => p.Mz).ToArray();
            var spectrum = new MzSpectrum(peaks.Select(p => p.Mz).ToArray(), peaks.Select(p => p.Height).ToArray(), false);
            return new MsDataScan(spectrum, scan + 1, 2, true, Polarity.Positive, scan, new MzRange(150, 1500), "synthetic", MZAnalyzerType.Orbitrap,
                peaks.Sum(p => p.Height), 10, null, $"scan={scan + 1}", isolationMZ: 500, isolationWidth: 25, dissociationType: DissociationType.HCD);
        }).ToArray();

    private static readonly (double Mz, double Intensity)[] SharedIons = [(300, 1), (400, 1), (500, 1)];
    private static readonly double[] SharedMzs = [300, 400, 500];

    private static (double Mz, Func<int, double> Height)[] Shared() => [(300, s => Peak(s)), (400, s => 0.8 * Peak(s)), (500, s => 0.6 * Peak(s))];

    private static double OwnIons(MsDataScan[] window, params (double Mz, double Intensity)[] own) =>
        DiaInterferenceRemoval.OwnIonCoElution(SharedIons.Concat(own).ToArray(), SharedMzs.Concat(own.Select(i => i.Mz)).ToArray(),
            SharedIons, SharedMzs, [window], PeakApexRt, Ppm, 3);

    [Test]
    public void OwnIonsOnTheSharedPeakCoElute()
    {
        var window = Window([.. Shared(), (900, s => 0.5 * Peak(s))]);

        Assert.That(OwnIons(window, (900, 1)), Is.GreaterThan(0.99));
    }

    [Test]
    public void OwnIonsThatAreNotThereDoNotCoElute()
    {
        Assert.That(OwnIons(Window(Shared()), (900, 1)), Is.EqualTo(0));
    }

    /// <summary>An own ion can be present and still belong to something else: its peak sits elsewhere in the same span.</summary>
    [Test]
    public void OwnIonsPeakingElsewhereDoNotCoElute()
    {
        var window = Window([.. Shared(), (900, s => Peak(s, apex: 3))]);

        Assert.That(OwnIons(window, (900, 1)), Is.LessThan(0.1));
    }

    /// <summary>A strong library ion that is missing says more than a weak one that is there.</summary>
    [Test]
    public void OwnIonsAreWeightedByTheLibrary()
    {
        var window = Window([.. Shared(), (900, s => 0.5 * Peak(s))]);

        Assert.That(OwnIons(window, (900, 1), (950, 3)), Is.EqualTo(0.25).Within(0.01));
    }

    /// <summary>With nothing of its own, a winner cannot be told apart from its rival at all.</summary>
    [Test]
    public void AWinnerWithNoOwnIonsDoesNotCoElute()
    {
        Assert.That(OwnIons(Window(Shared())), Is.EqualTo(0));
    }

    /// <summary>Shared means anywhere on the rival's whole list, scorable or not, because that signal cannot tell the two apart.</summary>
    [Test]
    public void AnIonAnywhereOnTheRivalsListIsNotTheWinnersOwn()
    {
        var window = Window([.. Shared(), (900, s => 0.5 * Peak(s))]);
        (double Mz, double Intensity)[] winner = [.. SharedIons, (900, 1)];

        double value = DiaInterferenceRemoval.OwnIonCoElution(winner, winner.Select(i => i.Mz).ToArray(), SharedIons, [.. SharedMzs, 900], [window], PeakApexRt, Ppm, 3);

        Assert.That(value, Is.EqualTo(0));
    }

    /// <summary>How much of a winner's library intensity is its own. Near 0 it has almost nothing to tell it from its rival.</summary>
    [Test]
    public void OwnIonShareIsTheLibraryIntensityNotOnTheRivalsList()
    {
        (double Mz, double Intensity)[] winner = [(300, 6), (400, 3), (900, 1)];

        Assert.That(DiaInterferenceRemoval.OwnIonShare(winner, [300, 400], Ppm), Is.EqualTo(0.1).Within(1e-9));
        Assert.That(DiaInterferenceRemoval.OwnIonShare(winner, [300, 400, 900], Ppm), Is.EqualTo(0));
        Assert.That(DiaInterferenceRemoval.OwnIonShare(winner, [], Ppm), Is.EqualTo(1));
        Assert.That(DiaInterferenceRemoval.OwnIonShare([], [300], Ppm), Is.EqualTo(0));
    }

    /// <summary>Shared signal found nowhere leaves nothing to measure against.</summary>
    [Test]
    public void NoSharedSignalIsNotANumber()
    {
        Assert.That(OwnIons(Window((900, s => Peak(s))), (900, 1)), Is.NaN);
    }

    /// <summary>A precursor can sit in overlapping windows. The one holding the shared signal is measured.</summary>
    [Test]
    public void TheWindowWithTheSharedSignalIsMeasured()
    {
        var empty = Window((950, s => Peak(s)));
        var real = Window([.. Shared(), (900, s => 0.5 * Peak(s))]);
        (double Mz, double Intensity)[] winner = [.. SharedIons, (900, 1)];

        double value = DiaInterferenceRemoval.OwnIonCoElution(winner, winner.Select(i => i.Mz).ToArray(), SharedIons, SharedMzs, [empty, real], PeakApexRt, Ppm, 3);

        Assert.That(value, Is.GreaterThan(0.99));
    }

    [Test]
    public void ArgumentsAreChecked()
    {
        Assert.Throws<ArgumentNullException>(() => DiaInterferenceRemoval.Remove(null!, _ => [], _ => 0, 1, 20, 4));
        Assert.That(DiaInterferenceRemoval.Remove([], _ => [], _ => 0, 1, 20, 4), Is.Empty);
    }
}
