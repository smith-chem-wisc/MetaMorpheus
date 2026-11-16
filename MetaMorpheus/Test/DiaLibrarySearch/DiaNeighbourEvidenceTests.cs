#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using EngineLayer.DiaLibrarySearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using StatisticalModels;

namespace Test.DiaLibrarySearch;

/// <summary>
/// How much of a precursor's fragment signal confident precursors eluting with it already explain, and what is left. On both benchmark runs decoys and false targets in a target's score
/// band share about twice as many ions with what elutes beside them as they do with precursors far away in retention time, while targets share no more than chance.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaNeighbourEvidenceTests
{
    private const double Ppm = 20;
    private const double RtTolerance = 2;
    private const int ApexHalfWidth = 3;
    private const int ApexScan = 10;
    private const string Sequence = "PEPTIDEK";

    // Four ions on one peak, in the library's proportions, so with nothing else the row reads a cosine of 1
    private static readonly (double Mz, double Intensity)[] Ions = [(300, 1.0), (400, 0.8), (500, 0.6), (600, 0.4)];
    private const double Sum = 2.8;

    private static double Peak(int scan, int apex = 10) => 100_000 * Math.Exp(-(scan - apex) * (scan - apex) / 8.0);

    private static MsDataScan[] Window(params (double Mz, Func<int, double> Height)[] ions) =>
        Enumerable.Range(0, 21).Select(scan =>
        {
            var peaks = ions.Select(i => (i.Mz, Height: i.Height(scan))).Where(p => p.Height > 0).OrderBy(p => p.Mz).ToArray();
            var spectrum = new MzSpectrum(peaks.Select(p => p.Mz).ToArray(), peaks.Select(p => p.Height).ToArray(), false);
            return new MsDataScan(spectrum, scan + 1, 2, true, Polarity.Positive, scan, new MzRange(150, 1500), "synthetic", MZAnalyzerType.Orbitrap,
                peaks.Sum(p => p.Height), 10, null, $"scan={scan + 1}", isolationMZ: 500, isolationWidth: 25, dissociationType: DissociationType.HCD);
        }).ToArray();

    private static (double Mz, Func<int, double> Height)[] OnThePeak() => Ions.Select(i => (i.Mz, (Func<int, double>)(s => i.Intensity * Peak(s)))).ToArray();

    // A neighbour holding the row's first ion, with three ions of its own at the given share of the peak
    private static DiaNeighbourEvidence.Neighbour Holding300(double apexRt = 10, string sequence = "NEIGHBOURK") =>
        new(apexRt, sequence, [1100, 1200, 1300, 300], [1, 1, 1, 1]);

    private static (double Mz, Func<int, double> Height)[] Its3Ions(double share) =>
        [(1100, s => share * Peak(s)), (1200, s => share * Peak(s)), (1300, s => share * Peak(s))];

    private static double[] Measure(MsDataScan[] window, DiaNeighbourEvidence.Neighbour[] neighbours, (double Mz, double Intensity)[]? ions = null, int apexScan = ApexScan, string sequence = Sequence)
    {
        ions ??= Ions;
        return DiaNeighbourEvidence.Measure(ions, ions.Select(i => i.Mz).ToArray(), sequence, neighbours, window, apexScan, RtTolerance, Ppm, ApexHalfWidth);
    }

    private static double Column(double[] values, string name) => values[Array.IndexOf(DiaNeighbourEvidence.ColumnNames, name)];

    [Test]
    public void WithNoNeighbourNothingIsExplained()
    {
        var read = Measure(Window(OnThePeak()), []);

        Assert.That(read, Has.Length.EqualTo(DiaNeighbourEvidence.ColumnNames.Length));
        Assert.That(Column(read, "ExplainedShare"), Is.EqualTo(0));
        Assert.That(Column(read, "ExplainedShareTop6"), Is.EqualTo(0));
        Assert.That(Column(read, "ResidualCosine"), Is.EqualTo(1).Within(1e-9));
        Assert.That(Column(read, "CleanIons"), Is.EqualTo(4));
        Assert.That(Column(read, "CleanLibraryShare"), Is.EqualTo(1));
        Assert.That(Column(read, "CleanCoElution"), Is.GreaterThan(0.99));
        Assert.That(Column(read, "SharedIons"), Is.EqualTo(0));
        Assert.That(Column(read, "Neighbours"), Is.EqualTo(0));
    }

    [Test]
    public void ANeighbourExplainsTheIonItHolds()
    {
        var read = Measure(Window([.. OnThePeak(), .. Its3Ions(1)]), [Holding300()]);

        Assert.That(Column(read, "ExplainedShare"), Is.EqualTo(1 / Sum).Within(1e-9));
        Assert.That(Column(read, "SharedIons"), Is.EqualTo(1));
        Assert.That(Column(read, "Neighbours"), Is.EqualTo(1));
        Assert.That(Column(read, "CleanIons"), Is.EqualTo(3));
        Assert.That(Column(read, "CleanLibraryShare"), Is.EqualTo(1.8 / Sum).Within(1e-9));
    }

    /// <summary>What the neighbour's own ions say it weighs here sets how much of the shared ion is its, however strong that ion reads.</summary>
    [Test]
    public void TheNeighbourIsScaledByTheIonsItDoesNotShare()
    {
        var read = Measure(Window([.. OnThePeak(), .. Its3Ions(0.25)]), [Holding300()]);

        Assert.That(Column(read, "ExplainedShare"), Is.EqualTo(0.25 / Sum).Within(1e-9));
        Assert.That(Column(read, "CleanIons"), Is.EqualTo(4), "a quarter explained still leaves the ion mostly its own");
    }

    /// <summary>The median, so one interfered ion of the neighbour's does not inflate it.</summary>
    [Test]
    public void OneInterferedIonDoesNotInflateTheNeighbour()
    {
        var window = Window([.. OnThePeak(), (1100, s => 9 * Peak(s)), (1200, s => 0.25 * Peak(s)), (1300, s => 0.25 * Peak(s))]);

        Assert.That(Column(Measure(window, [Holding300()]), "ExplainedShare"), Is.EqualTo(0.25 / Sum).Within(1e-9));
    }

    /// <summary>A near-isomer has too few ions of its own to weigh it by, so all its top ions do. Weighed by the two ions it has, a near-isomer would explain almost nothing.</summary>
    [Test]
    public void ANeighbourWithFewIonsOfItsOwnIsScaledOnAllItsTopIons()
    {
        DiaNeighbourEvidence.Neighbour nearIsomer = new(10, "NEIGHBOURK", [300, 400, 500, 1100, 1200], [1.0, 0.8, 0.6, 1.0, 1.0]);
        var window = Window([.. OnThePeak(), (1100, s => 0.01 * Peak(s)), (1200, s => 0.01 * Peak(s))]);

        Assert.That(Column(Measure(window, [nearIsomer]), "ExplainedShare"), Is.EqualTo(2.4 / Sum).Within(1e-9));
    }

    /// <summary>Three ions of its own are enough, even when they say it is not here.</summary>
    [Test]
    public void ThreeIonsOfItsOwnAreEnoughToScaleANeighbour()
    {
        DiaNeighbourEvidence.Neighbour neighbour = new(10, "NEIGHBOURK", [300, 400, 500, 1100, 1200, 1300], [1.0, 0.8, 0.6, 1.0, 1.0, 1.0]);

        Assert.That(Column(Measure(Window(OnThePeak()), [neighbour]), "ExplainedShare"), Is.EqualTo(0));
    }

    /// <summary>Only its six strongest ions weigh it.</summary>
    [Test]
    public void OnlyANeighboursTopSixIonsScaleIt()
    {
        DiaNeighbourEvidence.Neighbour neighbour = new(10, "NEIGHBOURK", [300, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000],
            Enumerable.Repeat(1.0, 11).ToArray());
        // Present among the weak five, absent among the strong: weighing all of them would give the middle one, on the peak
        var late = new[] { 1600.0, 1700, 1800, 1900, 2000 }.Select(mz => (mz, (Func<int, double>)(s => Peak(s)))).ToArray();
        var read = Measure(Window([.. OnThePeak(), .. late]), [neighbour]);

        Assert.That(Column(read, "ExplainedShare"), Is.EqualTo(0));
    }

    /// <summary>Its other charge states are the same peptide, so their signal is its own.</summary>
    [Test]
    public void ItsOwnChargeStatesExplainNothing()
    {
        var read = Measure(Window([.. OnThePeak(), .. Its3Ions(1)]), [Holding300(sequence: Sequence)]);

        Assert.That(Column(read, "ExplainedShare"), Is.EqualTo(0));
        Assert.That(Column(read, "Neighbours"), Is.EqualTo(0));
        Assert.That(Column(read, "SharedIons"), Is.EqualTo(0));
    }

    [TestCase(8, 1)]
    [TestCase(12, 1)]
    [TestCase(7.9, 0)]
    [TestCase(12.1, 0)]
    public void ANeighbourCountsWithinTheRtTolerance(double apexRt, int expected)
    {
        var read = Measure(Window([.. OnThePeak(), .. Its3Ions(1)]), [Holding300(apexRt)]);

        Assert.That(Column(read, "Neighbours"), Is.EqualTo(expected));
        Assert.That(Column(read, "SharedIons"), Is.EqualTo(expected));
    }

    /// <summary>Neighbours arrive sorted by apex, and every one in tolerance is seen, not only the first.</summary>
    [Test]
    public void EveryNeighbourInToleranceIsCounted()
    {
        var read = Measure(Window([.. OnThePeak(), .. Its3Ions(1)]), [Holding300(1), Holding300(8), Holding300(10), Holding300(12), Holding300(18)]);

        Assert.That(Column(read, "Neighbours"), Is.EqualTo(3));
    }

    /// <summary>Holding the ion is not enough: a neighbour whose own ions are not there explains nothing, though the sharing is counted.</summary>
    [Test]
    public void ANeighbourWhoseIonsAreAbsentExplainsNothing()
    {
        var read = Measure(Window(OnThePeak()), [Holding300()]);

        Assert.That(Column(read, "ExplainedShare"), Is.EqualTo(0));
        Assert.That(Column(read, "SharedIons"), Is.EqualTo(1));
        Assert.That(Column(read, "CleanIons"), Is.EqualTo(4));
    }

    [Test]
    public void NeighboursCannotExplainMoreThanIsThere()
    {
        var read = Measure(Window([.. OnThePeak(), .. Its3Ions(5)]), [Holding300()]);

        Assert.That(Column(read, "ExplainedShare"), Is.EqualTo(1 / Sum).Within(1e-9));
    }

    [Test]
    public void NeighboursExplainTogether()
    {
        DiaNeighbourEvidence.Neighbour other = new(10, "OTHERK", [1400, 1500, 1600, 300], [1, 1, 1, 1]);
        var window = Window([.. OnThePeak(), .. Its3Ions(0.3), (1400, s => 0.3 * Peak(s)), (1500, s => 0.3 * Peak(s)), (1600, s => 0.3 * Peak(s))]);
        var read = Measure(window, [Holding300(), other]);

        Assert.That(Column(read, "ExplainedShare"), Is.EqualTo(0.6 / Sum).Within(1e-9));
        Assert.That(Column(read, "CleanIons"), Is.EqualTo(3), "together they explain more than half of that ion");
    }

    /// <summary>A neighbour's scale multiplies its library intensity at the shared ion, not 1.</summary>
    [Test]
    public void AWeakIonOfTheNeighbourExplainsLittle()
    {
        DiaNeighbourEvidence.Neighbour neighbour = new(10, "NEIGHBOURK", [1100, 1200, 1300, 300], [1, 1, 1, 0.1]);

        Assert.That(Column(Measure(Window([.. OnThePeak(), .. Its3Ions(1)]), [neighbour]), "ExplainedShare"), Is.EqualTo(0.1 / Sum).Within(1e-9));
    }

    /// <summary>The row's whole library list decides what the neighbour shares when it is weighed, so an ion the row does not score cannot weigh it either.</summary>
    [Test]
    public void AnIonOnTheRowsLibraryListDoesNotWeighTheNeighbour()
    {
        DiaNeighbourEvidence.Neighbour neighbour = new(10, "NEIGHBOURK", [1100, 1200, 1300, 1400, 1500, 300], [1, 1, 1, 1, 1, 1]);
        var window = Window([.. OnThePeak(), (1100, s => 0.1 * Peak(s)), (1200, s => 9 * Peak(s)), (1300, s => 0.1 * Peak(s)), (1400, s => 9 * Peak(s)), (1500, s => 0.25 * Peak(s))]);
        double[] allMzs = Ions.Select(i => i.Mz).Append(1200).Append(1400).ToArray();

        // 1200 and 1400 are the row's too, leaving 0.1, 0.1 and 0.25, with a median of 0.1. Read as free, they would make it 0.25
        var read = DiaNeighbourEvidence.Measure(Ions, allMzs, Sequence, [neighbour], window, ApexScan, RtTolerance, Ppm, ApexHalfWidth);

        Assert.That(Column(read, "ExplainedShare"), Is.EqualTo(0.1 / Sum).Within(1e-9));
    }

    [Test]
    public void TheResidualIsWhatTheNeighboursLeave()
    {
        var read = Measure(Window([.. OnThePeak(), .. Its3Ions(1)]), [Holding300()]);

        // Library 1, 0.8, 0.6, 0.4 against what is left, 0, 0.8, 0.6, 0.4, on the square-root scale
        Assert.That(Column(read, "ResidualCosine"), Is.EqualTo(1.8 / Math.Sqrt(Sum * 1.8)).Within(1e-9));
    }

    /// <summary>An ion on another peak lowers the co-elution until a neighbour takes it, and then it is not read.</summary>
    [Test]
    public void CleanCoElutionReadsOnlyTheCleanIons()
    {
        var ions = new (double Mz, Func<int, double> Height)[] { (300, s => Peak(s, apex: 8)), (400, s => 0.8 * Peak(s)), (500, s => 0.6 * Peak(s)), (600, s => 0.4 * Peak(s)) };
        double all = Column(Measure(Window(ions), []), "CleanCoElution");
        var window = Window([.. ions, (1100, s => Peak(s, apex: 8)), (1200, s => Peak(s, apex: 8)), (1300, s => Peak(s, apex: 8))]);
        var read = Measure(window, [Holding300()]);

        Assert.That(Column(read, "CleanIons"), Is.EqualTo(3));
        Assert.That(Column(read, "CleanCoElution"), Is.GreaterThan(0.99));
        Assert.That(all, Is.LessThan(0.9));
    }

    /// <summary>Each ion is compared with the others, not with a sum holding itself, or an ion alone on its peak would read well.</summary>
    [Test]
    public void AnIonIsNotComparedWithItself()
    {
        var read = Measure(Window((300, s => Peak(s)), (400, s => 0.1 * Peak(s, apex: 7))), [], ions: [(300, 1), (400, 1)]);

        Assert.That(Column(read, "CleanCoElution"), Is.LessThan(0.3));
    }

    [Test]
    public void CleanCoElutionIsWeightedByTheLibrary()
    {
        var ions = new (double Mz, Func<int, double> Height)[] { (300, s => Peak(s)), (400, s => Peak(s)), (500, s => Peak(s, apex: 6)) };
        double weakOff = Column(Measure(Window(ions), [], ions: [(300, 1), (400, 1), (500, 0.01)]), "CleanCoElution");
        double strongOff = Column(Measure(Window(ions), [], ions: [(300, 1), (400, 1), (500, 1)]), "CleanCoElution");

        Assert.That(weakOff, Is.GreaterThan(strongOff));
    }

    [Test]
    public void WithFewerThanTwoCleanIonsCleanCoElutionIsZero()
    {
        DiaNeighbourEvidence.Neighbour takesThree = new(10, "NEIGHBOURK", [1100, 1200, 1300, 300, 400, 500], [1, 1, 1, 1, 0.8, 0.6]);
        var read = Measure(Window([.. OnThePeak(), .. Its3Ions(1)]), [takesThree]);

        Assert.That(Column(read, "CleanIons"), Is.EqualTo(1));
        Assert.That(Column(read, "CleanCoElution"), Is.EqualTo(0));
    }

    /// <summary>An ion the row does not see at all is not clean, or a row with nothing there would read as all clean.</summary>
    [Test]
    public void AnIonThatIsNotThereIsNotClean()
    {
        var read = Measure(Window((300, s => Peak(s)), (400, s => 0.8 * Peak(s))), []);

        Assert.That(Column(read, "CleanIons"), Is.EqualTo(2));
        Assert.That(Column(read, "CleanLibraryShare"), Is.EqualTo(1.8 / Sum).Within(1e-9));
    }

    /// <summary>Half explained is no longer its own.</summary>
    [Test]
    public void AnIonExactlyHalfExplainedIsNotClean()
    {
        var read = Measure(Window([.. OnThePeak(), .. Its3Ions(0.5)]), [Holding300()]);

        Assert.That(Column(read, "CleanIons"), Is.EqualTo(3));
    }

    [Test]
    public void AnIonUnderHalfExplainedIsClean()
    {
        var read = Measure(Window([.. OnThePeak(), .. Its3Ions(0.45)]), [Holding300()]);

        Assert.That(Column(read, "CleanIons"), Is.EqualTo(4));
    }

    /// <summary>The trees cannot be handed a value that is not a number.</summary>
    [Test]
    public void ARowWithNoSignalReadsAsNothing()
    {
        var read = Measure(Window((1100, s => Peak(s))), [Holding300()]);

        Assert.That(read, Has.All.Matches<double>(double.IsFinite));
        Assert.That(Column(read, "ExplainedShare"), Is.EqualTo(0));
        Assert.That(Column(read, "ExplainedShareTop6"), Is.EqualTo(0));
        Assert.That(Column(read, "ResidualCosine"), Is.EqualTo(0));
        Assert.That(Column(read, "CleanIons"), Is.EqualTo(0));
        Assert.That(Column(read, "CleanLibraryShare"), Is.EqualTo(0));
        Assert.That(Column(read, "CleanCoElution"), Is.EqualTo(0));
    }

    [Test]
    public void ARowWithNoIonsReadsAsNothing()
    {
        var read = Measure(Window(OnThePeak()), [Holding300()], ions: []);

        Assert.That(read, Has.All.Matches<double>(double.IsFinite));
        Assert.That(Column(read, "CleanLibraryShare"), Is.EqualTo(0));
    }

    [Test]
    public void TheTopSixShareReadsTheStrongestSix()
    {
        (double Mz, double Intensity)[] ions = [(300, 1), (400, 1), (500, 1), (600, 1), (650, 1), (700, 1), (800, 1), (900, 1)];
        var peaks = ions.Select(i => (i.Mz, (Func<int, double>)(s => Peak(s)))).ToArray();
        DiaNeighbourEvidence.Neighbour holdsTheSeventh = new(10, "NEIGHBOURK", [1100, 1200, 1300, 800], [1, 1, 1, 1]);
        var read = Measure(Window([.. peaks, .. Its3Ions(1)]), [holdsTheSeventh], ions: ions);

        Assert.That(Column(read, "ExplainedShare"), Is.EqualTo(1.0 / 8).Within(1e-9));
        Assert.That(Column(read, "ExplainedShareTop6"), Is.EqualTo(0));

        DiaNeighbourEvidence.Neighbour holdsTheSixth = new(10, "NEIGHBOURK", [1100, 1200, 1300, 700], [1, 1, 1, 1]);
        Assert.That(Column(Measure(Window([.. peaks, .. Its3Ions(1)]), [holdsTheSixth], ions: ions), "ExplainedShareTop6"), Is.EqualTo(1.0 / 6).Within(1e-9));
    }

    [TestCase(0)]
    [TestCase(20)]
    public void AnApexAtTheEdgeOfTheRunIsRead(int apexScan)
    {
        var flat = Ions.Select(i => (i.Mz, (Func<int, double>)(s => i.Intensity * 1000))).ToArray();
        var read = Measure(Window(flat), [], apexScan: apexScan);

        Assert.That(Column(read, "ResidualCosine"), Is.EqualTo(1).Within(1e-9));
        Assert.That(Column(read, "CleanIons"), Is.EqualTo(4));
    }

    /// <summary>The apex is read over its three scans: a neighbour's ions that dip on the apex scan still weigh it.</summary>
    [Test]
    public void TheApexIsReadOverThreeScans()
    {
        static Func<int, double> Beside(double height) => s => s is 9 or 11 ? height : 0;
        var window = Window((300, s => s == 10 ? 100 : 0), (1100, Beside(50)), (1200, Beside(50)), (1300, Beside(50)));
        var read = Measure(window, [Holding300()], ions: [(300, 1), (2000, 1)]);

        Assert.That(Column(read, "ExplainedShare"), Is.EqualTo(1));
    }

    /// <summary>Peaks two scans off the apex are outside the three.</summary>
    [Test]
    public void ThreeScansAreAllThatIsRead()
    {
        static Func<int, double> TwoOff(double height) => s => s is 8 or 12 ? height : 0;
        var window = Window((300, s => s == 10 ? 100 : 0), (1100, TwoOff(50)), (1200, TwoOff(50)), (1300, TwoOff(50)));

        Assert.That(Column(Measure(window, [Holding300()], ions: [(300, 1), (2000, 1)]), "ExplainedShare"), Is.EqualTo(0));
    }

    [Test]
    public void TheToleranceDecidesWhatIsShared()
    {
        DiaNeighbourEvidence.Neighbour inside = new(10, "NEIGHBOURK", [1100, 1200, 1300, 300.005], [1, 1, 1, 1]);
        DiaNeighbourEvidence.Neighbour outside = new(10, "NEIGHBOURK", [1100, 1200, 1300, 300.007], [1, 1, 1, 1]);
        var window = Window([.. OnThePeak(), .. Its3Ions(1)]);

        Assert.That(Column(Measure(window, [inside]), "SharedIons"), Is.EqualTo(1), "16.7 ppm");
        Assert.That(Column(Measure(window, [outside]), "SharedIons"), Is.EqualTo(0), "23.3 ppm");
    }

    // ---- The stack ----

    private static (double[] Scores, double[][] Columns, bool[] IsDecoy, string[] Keys) Rows(int count = 2000, Func<int, Random, bool, double>? column = null, int decoyEvery = 2)
    {
        var random = new Random(7);
        column ??= (_, r, decoy) => (decoy ? 0 : 1) + 0.05 * r.NextDouble();
        var isDecoy = Enumerable.Range(0, count).Select(i => i % decoyEvery == 0).ToArray();
        return (Enumerable.Range(0, count).Select(_ => random.NextDouble()).ToArray(), isDecoy.Select((d, i) => new[] { column(i, random, d) }).ToArray(), isDecoy,
            Enumerable.Range(0, count).Select(i => $"PEPTIDE{i}K").ToArray());
    }

    [Test]
    public void TheStackKeepsTheScoreValues()
    {
        var (scores, columns, isDecoy, keys) = Rows();
        double[] stacked = DiaNeighbourEvidence.Stack(scores, columns, isDecoy, keys, 0);

        Assert.That(stacked.Order(), Is.EqualTo(scores.Order()));
    }

    [Test]
    public void TheStackOrdersByTheEvidence()
    {
        var (scores, columns, isDecoy, keys) = Rows();
        double[] stacked = DiaNeighbourEvidence.Stack(scores, columns, isDecoy, keys, 0);
        int[] top = Enumerable.Range(0, scores.Length).OrderByDescending(i => stacked[i]).Take(500).ToArray();

        Assert.That(top.Count(i => !isDecoy[i]), Is.GreaterThan(475), "the column separates them, the score does not");
    }

    /// <summary>The score goes in as a column too, so evidence that says nothing leaves the order mostly as it was.</summary>
    [Test]
    public void TheStackLearnsFromTheScore()
    {
        var random = new Random(3);
        var isDecoy = Enumerable.Range(0, 2000).Select(i => i % 2 == 0).ToArray();
        double[] scores = isDecoy.Select(d => (d ? 0 : 1) + 0.05 * random.NextDouble()).ToArray();
        double[][] noise = isDecoy.Select(_ => new[] { random.NextDouble() }).ToArray();
        double[] stacked = DiaNeighbourEvidence.Stack(scores, noise, isDecoy, Enumerable.Range(0, 2000).Select(i => $"PEPTIDE{i}K").ToArray(), 0);
        int[] top = Enumerable.Range(0, 2000).OrderByDescending(i => stacked[i]).Take(500).ToArray();

        Assert.That(top.Count(i => !isDecoy[i]), Is.GreaterThan(475));
    }

    /// <summary>A column naming each row's group, whose label goes with it, is perfect in-sample. Out of fold the group is unseen, so it can only mislead.</summary>
    [Test]
    public void EachRowIsScoredByTreesThatNeverSawIt()
    {
        var random = new Random(7);
        int[] group = Enumerable.Range(0, 2000).Select(i => i % 20).ToArray();
        bool[] isDecoy = group.Select(g => g % 2 == 0).ToArray();
        double[] scores = group.Select(_ => random.NextDouble()).ToArray();
        double[] stacked = DiaNeighbourEvidence.Stack(scores, group.Select(g => new double[] { g }).ToArray(), isDecoy, group.Select(g => $"PEPTIDE{g:00}K").ToArray(), 0);
        int[] top = Enumerable.Range(0, 2000).OrderByDescending(i => stacked[i]).Take(500).ToArray();

        Assert.That(top.Count(i => !isDecoy[i]), Is.LessThan(375));
    }

    /// <summary>Trees are coarse at the very top: on the held-out run, unguarded, a few decoys climbing there cost 28% of the targets at q 0.1% while nothing else moved.</summary>
    [Test]
    public void TheVeryTopKeepsItsPlace()
    {
        var (scores, columns, isDecoy, keys) = Rows();
        double keepAtOrAbove = scores.OrderDescending().ElementAt(200);
        double[] stacked = DiaNeighbourEvidence.Stack(scores, columns, isDecoy, keys, 0, keepAtOrAbove: keepAtOrAbove);
        int[] rest = Enumerable.Range(0, scores.Length).Where(i => scores[i] < keepAtOrAbove).ToArray();

        Assert.That(Enumerable.Range(0, scores.Length).Where(i => scores[i] >= keepAtOrAbove).Select(i => stacked[i]),
            Is.EqualTo(Enumerable.Range(0, scores.Length).Where(i => scores[i] >= keepAtOrAbove).Select(i => scores[i])));
        Assert.That(rest.Select(i => stacked[i]).Max(), Is.LessThan(keepAtOrAbove), "nothing climbs past them");
        Assert.That(rest.Count(i => stacked[i] != scores[i]), Is.GreaterThan(rest.Length / 2), "the rest are still re-ranked");
        int[] topOfRest = rest.OrderByDescending(i => stacked[i]).Take(400).ToArray();
        Assert.That(topOfRest.Count(i => !isDecoy[i]), Is.GreaterThan(380), "by the evidence");
    }

    [Test]
    public void WithTooFewDecoysTheStackLeavesTheScoresAlone()
    {
        var (scores, columns, _, keys) = Rows();
        bool[] fewDecoys = scores.Select((_, i) => i % 200 == 0).ToArray();

        Assert.That(DiaNeighbourEvidence.Stack(scores, columns, fewDecoys, keys, 0), Is.EqualTo(scores));
    }

    [Test]
    public void TheStackRefusesListsOfDifferentLengths()
    {
        var (scores, columns, isDecoy, keys) = Rows();

        Assert.That(() => DiaNeighbourEvidence.Stack(scores, columns.Skip(1).ToArray(), isDecoy, keys, 0), Throws.ArgumentException);
        Assert.That(() => DiaNeighbourEvidence.Stack(scores, columns, isDecoy.Skip(1).ToArray(), keys, 0), Throws.ArgumentException);
        Assert.That(() => DiaNeighbourEvidence.Stack(scores, columns, isDecoy, keys.Skip(1).ToArray(), 0), Throws.ArgumentException);
    }

    // ---- The whole step ----

    // 300 confident targets over two isolation windows and 21 scans, then the row under test, then decoys below all of them
    private static DiaPrecursorMatch Match(int index, double score, bool isDecoy = false, double rt = 10, double mz = 500) =>
        new(index, $"PEPTIDE{index}K", 2, mz, isDecoy, new Irt(0), new RtMinutes(rt), new Irt(0), score, [score]);

    private static List<DiaPrecursorMatch> Crowd() =>
        Enumerable.Range(1000, 300).Select(i => Match(i, 10 + i / 1000.0, rt: 1)).ToList();

    private static List<DiaPrecursorMatch> Tail() =>
        Enumerable.Range(2000, 200).Select(i => Match(i, -10 - i / 1000.0, isDecoy: i % 2 == 0, rt: 1)).ToList();

    private static DiaNeighbourEvidence.PrecursorIons IonsOf(int index) =>
        new(Enumerable.Range(0, 4).Select(i => (5000.0 + index + i, 1.0)).ToList(), [], $"BASE{index}");

    [Test]
    public void TheColumnsNeverSeeADecoyOrAnUnconfidentTargetAsANeighbour()
    {
        // 20 near-top decoys bring the q of everything under them over 1% and leave it under 10%, so the two rows at the end are measured but not confident
        var matches = Crowd();
        matches.AddRange(Enumerable.Range(3000, 20).Select(i => Match(i, 5 + i / 10000.0, isDecoy: true, rt: 10)));
        var unconfident = Match(1, 4, rt: 10);
        var row = Match(0, 3, rt: 10);
        matches.Add(unconfident);
        matches.Add(row);
        matches.AddRange(Tail());
        var confidentRival = Match(2, 10, rt: 10);
        matches.Add(confidentRival);

        var window = Window([.. OnThePeak(), .. Its3Ions(1)]);
        DiaNeighbourEvidence.PrecursorIons Of(int index) => index switch
        {
            0 => new(Ions.ToList(), Ions.Select(i => i.Mz).ToList(), "ROWK"),
            1 or >= 3000 => new([(1100, 1), (1200, 1), (1300, 1), (300, 1)], [1100, 1200, 1300, 300], $"RIVAL{index}K"),
            2 => new([(1100, 1), (1200, 1), (1300, 1), (400, 1)], [1100, 1200, 1300, 400], "CONFIDENTK"),
            _ => IonsOf(index),
        };

        var columns = DiaNeighbourEvidence.MeasureAll(matches, Of, [(490, 510)], [window], RtTolerance, Ppm, ApexHalfWidth, true, 12);
        double[] read = columns[matches.IndexOf(row)]!;

        Assert.That(Column(read, "Neighbours"), Is.EqualTo(1), "only the confident target");
        Assert.That(Column(read, "SharedIons"), Is.EqualTo(1));
        Assert.That(Column(read, "ExplainedShare"), Is.EqualTo(0.8 / Sum).Within(1e-9));
    }

    [Test]
    public void ARowIsMeasuredWithItsTopIonsAndNeighboursWithAllTheirs()
    {
        var row = Match(0, 11, rt: 10);
        var neighbour = Match(1, 12, rt: 10);
        var matches = Crowd();
        matches.Add(row);
        matches.Add(neighbour);
        matches.AddRange(Tail());
        var window = Window([.. OnThePeak(), .. Its3Ions(1), (700, s => Peak(s))]);
        DiaNeighbourEvidence.PrecursorIons Of(int index) => index switch
        {
            // Its fifth ion is past the top four, so the neighbour holding it shares nothing that is read
            0 => new([.. Ions, (700, 0.1)], [300, 400, 500, 600, 700], "ROWK"),
            1 => new([(1100, 1), (1200, 1), (1300, 1), (600, 0.1), (700, 0.05)], [1100, 1200, 1300, 600, 700], "NEIGHBOURK"),
            _ => IonsOf(index),
        };

        double[] read = DiaNeighbourEvidence.MeasureAll(matches, Of, [(490, 510)], [window], RtTolerance, Ppm, ApexHalfWidth, true, 4)[matches.IndexOf(row)]!;

        Assert.That(Column(read, "SharedIons"), Is.EqualTo(1), "the neighbour's fourth ion is one of the row's top four; its fifth is on the row's list but past them");
        Assert.That(Column(read, "ExplainedShare"), Is.EqualTo(0.1 / Sum).Within(1e-9));
    }

    /// <summary>A precursor in two windows is read in the window holding its apex scan, and its neighbours are that window's.</summary>
    [Test]
    public void ARowIsReadInTheWindowHoldingItsApex()
    {
        var row = Match(0, 11, rt: 10.0, mz: 500);
        var neighbour = Match(1, 12, rt: 10.0, mz: 520);
        var matches = Crowd();
        matches.Add(row);
        matches.Add(neighbour);
        matches.AddRange(Tail());
        var holding = Window([.. OnThePeak(), .. Its3Ions(1)]);
        // The other window's scans fall between this one's, and hold nothing
        var other = Enumerable.Range(0, 21).Select(scan => new MsDataScan(new MzSpectrum(Array.Empty<double>(), Array.Empty<double>(), false), scan + 101, 2, true, Polarity.Positive, scan + 0.5, new MzRange(150, 1500), "synthetic", MZAnalyzerType.Orbitrap,
            0, 10, null, $"scan={scan + 101}", isolationMZ: 500, isolationWidth: 25, dissociationType: DissociationType.HCD)).ToArray();
        DiaNeighbourEvidence.PrecursorIons Of(int index) => index switch
        {
            0 => new(Ions.ToList(), Ions.Select(i => i.Mz).ToList(), "ROWK"),
            1 => new([(1100, 1), (1200, 1), (1300, 1), (300, 1)], [1100, 1200, 1300, 300], "NEIGHBOURK"),
            _ => IonsOf(index),
        };

        // Listed second so that taking the first window holding the m/z would read the empty one; the neighbour's m/z only falls in the holding window
        var columns = DiaNeighbourEvidence.MeasureAll(matches, Of, [(490, 510), (490, 530)], [other, holding], RtTolerance, Ppm, ApexHalfWidth, true, 12);
        double[] read = columns[matches.IndexOf(row)]!;

        Assert.That(Column(read, "ResidualCosine"), Is.GreaterThan(0));
        Assert.That(Column(read, "Neighbours"), Is.EqualTo(1));
    }

    /// <summary>A precursor isolated elsewhere puts nothing in this window's scans, however well it would fit. A window holds its low bound and not its high one.</summary>
    [TestCase(490.0, 1)]
    [TestCase(509.99, 1)]
    [TestCase(510.0, 0)]
    [TestCase(489.99, 0)]
    public void ANeighbourCountsOnlyInItsOwnWindow(double neighbourMz, int expected)
    {
        var row = Match(0, 11, rt: 10, mz: 500);
        var matches = Crowd();
        matches.Add(row);
        matches.Add(Match(1, 12, rt: 10, mz: neighbourMz));
        matches.AddRange(Tail());
        DiaNeighbourEvidence.PrecursorIons Of(int index) => index switch
        {
            0 => new(Ions.ToList(), Ions.Select(i => i.Mz).ToList(), "ROWK"),
            1 => new([(1100, 1), (1200, 1), (1300, 1), (300, 1)], [1100, 1200, 1300, 300], "NEIGHBOURK"),
            _ => IonsOf(index),
        };

        double[] read = DiaNeighbourEvidence.MeasureAll(matches, Of, [(490, 510)], [Window([.. OnThePeak(), .. Its3Ions(1)])], RtTolerance, Ppm, ApexHalfWidth, true, 12)[matches.IndexOf(row)]!;

        Assert.That(Column(read, "Neighbours"), Is.EqualTo(expected));
    }

    [Test]
    public void OnlyRowsDownToTheProvisionalFloorAreMeasured()
    {
        var matches = Crowd();
        matches.AddRange(Tail());
        var window = Window(OnThePeak());

        var columns = DiaNeighbourEvidence.MeasureAll(matches, IonsOf, [(490, 510)], [window], RtTolerance, Ppm, ApexHalfWidth, true, 12);

        Assert.That(columns[0], Is.Not.Null);
        Assert.That(columns[^1], Is.Null);
    }

    [Test]
    public void ApplyMovesOnlyRowsAboveTheFloor()
    {
        var random = new Random(5);
        // One decoy in 20 keeps every one of these under the provisional q; the decoys after them pass it
        var matches = Enumerable.Range(0, 2000).Select(i => Match(i, random.NextDouble() + 10, isDecoy: i % 20 == 0, rt: 1 + i % 19)).ToList();
        matches.AddRange(Tail().Select(m => m with { IsDecoy = true }));
        var window = Window(OnThePeak());
        DiaNeighbourEvidence.PrecursorIons Of(int index) => new(Ions.ToList(), Ions.Select(i => i.Mz).ToList(), $"BASE{index}");
        double[] q = TargetDecoyQValues.Compute(matches.Select(m => m.Score).ToArray(), matches.Select(m => m.IsDecoy).ToArray());
        double floor = matches.Where((_, i) => q[i] <= DiaInterferenceRemoval.ProvisionalQValue).Min(m => m.Score);

        var lines = new List<string>();
        var applied = DiaNeighbourEvidence.Apply(matches, Of, [(490, 510)], [window], RtTolerance, Ppm, ApexHalfWidth, true, 12, 0, 1, lines.Add);

        Assert.That(applied.Select(m => m.PrecursorIndex), Is.EqualTo(matches.Select(m => m.PrecursorIndex)), "same rows, same order");
        var below = Enumerable.Range(0, matches.Count).Where(i => matches[i].Score < floor).ToList();
        var above = Enumerable.Range(0, matches.Count).Where(i => matches[i].Score >= floor).ToList();
        Assert.That(below, Is.Not.Empty, "the case needs rows under the floor");
        Assert.That(below.Select(i => applied[i].Score), Is.EqualTo(below.Select(i => matches[i].Score)));
        Assert.That(above.Select(i => applied[i].Score).Order(), Is.EqualTo(above.Select(i => matches[i].Score).Order()), "the values stay the classifier's own");
        Assert.That(above.Count(i => applied[i].Score != matches[i].Score), Is.GreaterThan(above.Count / 2), "the order must really change");
        Assert.That(lines, Has.One.Contains("neighbour evidence"));
    }

    /// <summary>Rows at or under a provisional q of 0.1% keep their scores, and the trees still learn from them.</summary>
    [Test]
    public void ApplyKeepsRowsAtTheVeryTopInPlace()
    {
        var random = new Random(5);
        // 1,500 targets alone at the top reach q 0.1%; below them one decoy in 20, down to the floor
        var matches = Enumerable.Range(0, 1500).Select(i => Match(i + 5000, 100 + i / 100.0, rt: 1 + i % 19)).ToList();
        matches.AddRange(Enumerable.Range(0, 2000).Select(i => Match(i, random.NextDouble() + 10, isDecoy: i % 20 == 0, rt: 1 + i % 19)));
        matches.AddRange(Tail().Select(m => m with { IsDecoy = true }));
        DiaNeighbourEvidence.PrecursorIons Of(int index) => new(Ions.ToList(), Ions.Select(i => i.Mz).ToList(), $"BASE{index}");
        double[] q = TargetDecoyQValues.Compute(matches.Select(m => m.Score).ToArray(), matches.Select(m => m.IsDecoy).ToArray());
        var kept = Enumerable.Range(0, matches.Count).Where(i => q[i] <= DiaNeighbourEvidence.KeptQValue).ToList();
        Assert.That(kept, Has.Count.GreaterThan(1000), "the case needs rows that q");
        Assert.That(q[1500 + 1], Is.GreaterThan(DiaNeighbourEvidence.KeptQValue), "the case needs rows past it");

        var applied = DiaNeighbourEvidence.Apply(matches, Of, [(490, 510)], [Window(OnThePeak())], RtTolerance, Ppm, ApexHalfWidth, true, 12, 0);

        Assert.That(kept.Select(i => applied[i].Score), Is.EqualTo(kept.Select(i => matches[i].Score)));
        Assert.That(Enumerable.Range(1500, 2000).Count(i => applied[i].Score != matches[i].Score), Is.GreaterThan(1000), "the rest are still re-ranked");
    }

    [Test]
    public void ApplyWithNothingAboveTheFloorChangesNothing()
    {
        var matches = Enumerable.Range(0, 10).Select(i => Match(i, -i, isDecoy: true)).ToList();

        var applied = DiaNeighbourEvidence.Apply(matches, IonsOf, [(490, 510)], [Window(OnThePeak())], RtTolerance, Ppm, ApexHalfWidth, true, 12, 0);

        Assert.That(applied.Select(m => m.Score), Is.EqualTo(matches.Select(m => m.Score)));
    }

    /// <summary>A precursor no window holds, or one the library has lost, has nothing to read and stays in the stack on zeros.</summary>
    [Test]
    public void ARowWithNoWindowOrNoIonsReadsAsZeros()
    {
        var matches = Crowd();
        matches.Add(Match(0, 11, mz: 900));
        matches.Add(Match(1, 11.5, rt: 10));
        matches.AddRange(Tail());
        DiaNeighbourEvidence.PrecursorIons? Of(int index) => index == 1 ? null : IonsOf(index);

        var columns = DiaNeighbourEvidence.MeasureAll(matches, Of, [(490, 510)], [Window(OnThePeak())], RtTolerance, Ppm, ApexHalfWidth, true, 12);

        Assert.That(columns[matches.FindIndex(m => m.PrecursorIndex == 0)], Is.All.EqualTo(0));
        Assert.That(columns[matches.FindIndex(m => m.PrecursorIndex == 1)], Is.All.EqualTo(0));
    }
}
