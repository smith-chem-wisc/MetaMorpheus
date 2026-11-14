#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using EngineLayer.DiaLibrarySearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;

namespace Test.DiaLibrarySearch;

/// <summary>
/// What a run's confident targets say of the library's predicted intensities. On both benchmark runs y3 reads about twice what the library predicts, falling steadily to about 0.6x by y14, and a
/// table fitted on half of the targets raises the other half's 10th-percentile cosine by 0.03 and 0.06.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaFragmentIntensityModelTests
{
    internal static MslFragmentIon Ion(ProductType type, int number, float intensity = 1, float mz = 500, int charge = 1, double neutralLoss = 0) =>
        new() { ProductType = type, FragmentNumber = number, Intensity = intensity, Mz = mz, Charge = charge, NeutralLoss = neutralLoss };

    /// <summary>One anchor per item: ions with their library intensities, each read at the given multiple of it.</summary>
    private static DiaFragmentIntensityModel Fit(int anchors, params (MslFragmentIon Ion, double ObservedOverLibrary)[] ions) =>
        DiaFragmentIntensityModel.Fit(Enumerable.Range(0, anchors).Select(_ =>
            (IReadOnlyList<(MslFragmentIon Ion, double Observed)>)ions.Select(i => (i.Ion, i.Ion.Intensity * i.ObservedOverLibrary)).ToList()));

    internal static double Ratio(DiaFragmentIntensityModel model, MslFragmentIon a, MslFragmentIon b) => model.Multiplier(a) / model.Multiplier(b);

    [Test]
    public void AClassReadingStrongerThanTheLibraryPredictsIsRaised()
    {
        var y3 = Ion(ProductType.y, 3);
        var y6 = Ion(ProductType.y, 6);

        var model = Fit(40, (y3, 2), (y6, 1));

        Assert.That(Ratio(model, y3, y6), Is.EqualTo(2).Within(1e-6));
        Assert.That(model.Multiplier(y3), Is.GreaterThan(1));
        Assert.That(model.Multiplier(y6), Is.LessThan(1));
    }

    /// <summary>The ratio is of shares within each precursor, so a bright precursor counts as a faint one does.</summary>
    [Test]
    public void APrecursorsBrightnessDoesNotMatter()
    {
        var y3 = Ion(ProductType.y, 3);
        var y6 = Ion(ProductType.y, 6);
        var anchors = Enumerable.Range(0, 40).Select(i =>
            (IReadOnlyList<(MslFragmentIon Ion, double Observed)>)new List<(MslFragmentIon, double)> { (y3, 2.0 * Math.Pow(10, i % 5)), (y6, Math.Pow(10, i % 5)) });

        Assert.That(Ratio(DiaFragmentIntensityModel.Fit(anchors), y3, y6), Is.EqualTo(2).Within(1e-6));
    }

    [Test]
    public void TheLibraryIntensityIsTakenIntoAccount()
    {
        var y3 = Ion(ProductType.y, 3, intensity: 0.5f);
        var y6 = Ion(ProductType.y, 6, intensity: 1f);

        // Read in the library's proportions: nothing to correct
        Assert.That(Ratio(Fit(40, (y3, 1), (y6, 1)), y3, y6), Is.EqualTo(1).Within(1e-6));
    }

    [Test]
    public void ALargeSwingIsClamped()
    {
        var y3 = Ion(ProductType.y, 3);
        var y6 = Ion(ProductType.y, 6);

        var raised = Fit(40, (y3, 1000), (y6, 1));
        var lowered = Fit(40, (y3, 0.001), (y6, 1));

        Assert.That(raised.Multiplier(y3), Is.LessThanOrEqualTo(DiaFragmentIntensityModel.MaximumMultiplier));
        Assert.That(raised.Multiplier(y6), Is.GreaterThanOrEqualTo(1 / DiaFragmentIntensityModel.MaximumMultiplier));
        Assert.That(lowered.Multiplier(y3), Is.GreaterThanOrEqualTo(1 / DiaFragmentIntensityModel.MaximumMultiplier));
        Assert.That(DiaFragmentIntensityModel.MaximumMultiplier, Is.EqualTo(4));
    }

    [Test]
    public void AClassSeenTooRarelyStaysAsItIs()
    {
        var common = Ion(ProductType.y, 6);
        var alsoCommon = Ion(ProductType.y, 7);
        var rare = Ion(ProductType.y, 3);
        var anchors = Enumerable.Range(0, 40).Select(i => (IReadOnlyList<(MslFragmentIon Ion, double Observed)>)(i < DiaFragmentIntensityModel.MinimumIonsPerClass - 1
            ? new List<(MslFragmentIon, double)> { (common, 1.0), (alsoCommon, 1.0), (rare, 3.0) }
            : new List<(MslFragmentIon, double)> { (common, 1.0), (alsoCommon, 1.0) }));

        Assert.That(DiaFragmentIntensityModel.Fit(anchors).Multiplier(rare), Is.EqualTo(1));
    }

    [Test]
    public void AClassSeenOftenEnoughIsFitted()
    {
        var common = Ion(ProductType.y, 6);
        var rare = Ion(ProductType.y, 3);
        var anchors = Enumerable.Range(0, 40).Select(i => (IReadOnlyList<(MslFragmentIon Ion, double Observed)>)(i < DiaFragmentIntensityModel.MinimumIonsPerClass
            ? new List<(MslFragmentIon, double)> { (common, 1.0), (rare, 3.0) }
            : new List<(MslFragmentIon, double)> { (common, 1.0) }));

        Assert.That(DiaFragmentIntensityModel.Fit(anchors).Multiplier(rare), Is.Not.EqualTo(1));
    }

    /// <summary>A weak ion only counts when it happens to be seen, so its class would read high. Those ions do not fit the table, but the table still applies to them.</summary>
    [Test]
    public void WeakIonsDoNotFitTheTable()
    {
        var weak = Ion(ProductType.y, 3, intensity: 0.04f);
        var alsoWeak = Ion(ProductType.y, 6, intensity: 0.04f);
        var strong = Ion(ProductType.y, 7);
        var strongToo = Ion(ProductType.y, 8);

        var model = Fit(40, (weak, 50), (alsoWeak, 1), (strong, 1), (strongToo, 1));

        Assert.That(model.Multiplier(weak), Is.EqualTo(1));
        Assert.That(DiaFragmentIntensityModel.LibraryFloor, Is.EqualTo(0.05));
    }

    [Test]
    public void AnIonAtTheFloorFitsTheTable()
    {
        var atFloor = Ion(ProductType.y, 3, intensity: 0.05f);
        var strong = Ion(ProductType.y, 7);

        Assert.That(Fit(40, (atFloor, 3), (strong, 1)).Multiplier(atFloor), Is.Not.EqualTo(1));
    }

    /// <summary>Not seen says nothing of its ratio: the ion may be under the noise.</summary>
    [Test]
    public void AnIonNotSeenDoesNotCount()
    {
        var y3 = Ion(ProductType.y, 3);
        var y6 = Ion(ProductType.y, 6);
        var anchors = Enumerable.Range(0, 80).Select(i =>
            (IReadOnlyList<(MslFragmentIon Ion, double Observed)>)new List<(MslFragmentIon, double)> { (y3, i % 4 == 0 ? 0.0 : 2.0), (y6, 1.0) });

        Assert.That(Ratio(DiaFragmentIntensityModel.Fit(anchors), y3, y6), Is.EqualTo(2).Within(1e-6));
    }

    /// <summary>An ion missing from most apexes is still judged on the ones it is in. Counted as nothing, it would read as far too weak.</summary>
    [Test]
    public void AnIonMostlyNotSeenIsJudgedWhereItIsSeen()
    {
        var y3 = Ion(ProductType.y, 3);
        var y6 = Ion(ProductType.y, 6);
        var anchors = Enumerable.Range(0, 160).Select(i =>
            (IReadOnlyList<(MslFragmentIon Ion, double Observed)>)new List<(MslFragmentIon, double)> { (y3, i % 4 == 0 ? 2.0 : 0.0), (y6, 1.0) });

        Assert.That(DiaFragmentIntensityModel.Fit(anchors).Multiplier(y3), Is.GreaterThan(1));
    }

    [Test]
    public void BIonsAreInTheTable()
    {
        var b5 = Ion(ProductType.b, 5);
        var y6 = Ion(ProductType.y, 6);

        Assert.That(Fit(40, (b5, 3), (y6, 1)).Multiplier(b5), Is.GreaterThan(1));
    }

    [Test]
    public void TheFitsBoundsAreAsBenchmarked()
    {
        Assert.That(DiaFragmentIntensityModel.MinimumIonsPerClass, Is.EqualTo(30));
        Assert.That(DiaFragmentIntensityModel.LibraryFloor, Is.EqualTo(0.05));
        Assert.That(DiaFragmentIntensityModel.MaximumMultiplier, Is.EqualTo(4));
    }

    [Test]
    public void TheTableUsesTheMedian()
    {
        var y3 = Ion(ProductType.y, 3);
        var y6 = Ion(ProductType.y, 6);
        var anchors = Enumerable.Range(0, 101).Select(i =>
            (IReadOnlyList<(MslFragmentIon Ion, double Observed)>)new List<(MslFragmentIon, double)> { (y3, i < 10 ? 3.9 : 2.0), (y6, 1.0) });

        Assert.That(Ratio(DiaFragmentIntensityModel.Fit(anchors), y3, y6), Is.EqualTo(2).Within(1e-6));
    }

    [Test]
    public void ClassesDifferByTypeNumberAndCharge()
    {
        var y5 = Ion(ProductType.y, 5);
        var b5 = Ion(ProductType.b, 5);
        var y6 = Ion(ProductType.y, 6);
        var y5Twice = Ion(ProductType.y, 5, charge: 2);
        var reference = Ion(ProductType.y, 7);

        var model = Fit(40, (y5, 2), (b5, 0.5), (y6, 1.5), (y5Twice, 0.75), (reference, 1));

        double[] multipliers = [model.Multiplier(y5), model.Multiplier(b5), model.Multiplier(y6), model.Multiplier(y5Twice)];
        Assert.That(multipliers.Distinct().Count(), Is.EqualTo(4));
    }

    [Test]
    public void IonsFarFromTheTerminusShareAClass()
    {
        var y14 = Ion(ProductType.y, 14);
        var y20 = Ion(ProductType.y, 20);
        var y13 = Ion(ProductType.y, 13);
        var model = Fit(40, (y14, 0.5), (Ion(ProductType.y, 6), 1), (y13, 1));

        Assert.That(model.Multiplier(y20), Is.EqualTo(model.Multiplier(y14)));
        Assert.That(model.Multiplier(y13), Is.Not.EqualTo(model.Multiplier(y14)));
    }

    /// <summary>The fragment rule never scores these, and a neutral loss or an internal ion is its own kind of ion, so none of them fit the table or are moved by it.</summary>
    [TestCase(ProductType.y, 2, 500, 0)]
    [TestCase(ProductType.y, 5, 199, 0)]
    [TestCase(ProductType.y, 5, 500, -18)]
    [TestCase(ProductType.a, 5, 500, 0)]
    public void IonsOutsideTheTableStayAsTheyAre(ProductType type, int number, float mz, double neutralLoss)
    {
        var odd = Ion(type, number, mz: mz, neutralLoss: neutralLoss);
        var y6 = Ion(ProductType.y, 6);
        var y7 = Ion(ProductType.y, 7);

        Assert.That(Fit(40, (odd, 3), (y6, 1), (y7, 1)).Multiplier(odd), Is.EqualTo(1));
    }

    [Test]
    public void AnInternalIonStaysAsItIs()
    {
        var y5 = Ion(ProductType.y, 5);
        var inner = Ion(ProductType.y, 5);
        inner.SecondaryProductType = ProductType.b;
        inner.SecondaryFragmentNumber = 7;

        Assert.That(Fit(40, (y5, 3), (Ion(ProductType.y, 6), 1)).Multiplier(inner), Is.EqualTo(1));
    }

    [Test]
    public void ANeutralLossDoesNotFitTheClassOfItsIon()
    {
        var y5 = Ion(ProductType.y, 5);
        var y5Loss = Ion(ProductType.y, 5, neutralLoss: -18);
        var y6 = Ion(ProductType.y, 6);

        var model = Fit(40, (y5, 2), (y5Loss, 50), (y6, 1));

        Assert.That(Ratio(model, y5, y6), Is.EqualTo(2).Within(1e-6));
    }

    [Test]
    public void WithNoAnchorsNothingMoves()
    {
        var model = DiaFragmentIntensityModel.Fit([]);

        Assert.That(model.Multiplier(Ion(ProductType.y, 3)), Is.EqualTo(1));
    }

    [Test]
    public void AnAnchorWithOneIonSaysNothing()
    {
        var y3 = Ion(ProductType.y, 3);

        // Its share of itself is always 1, so it cannot say a class reads high
        Assert.That(Fit(40, (y3, 5)).Multiplier(y3), Is.EqualTo(1));
    }

    [Test]
    public void TheLibraryIntensityIsMultiplied()
    {
        var y3 = Ion(ProductType.y, 3, intensity: 0.5f);
        var y6 = Ion(ProductType.y, 6, intensity: 0.5f);
        var model = Fit(40, (y3, 2), (y6, 1));

        Assert.That(model.Intensity(y3), Is.EqualTo(0.5 * model.Multiplier(y3)).Within(1e-9));
    }

    // ---- Reading the anchors from the run ----

    private static MsDataScan[] Window(params (double Mz, Func<int, double> Height)[] ions) =>
        Enumerable.Range(0, 21).Select(scan =>
        {
            var peaks = ions.Select(i => (i.Mz, Height: i.Height(scan))).Where(p => p.Height > 0).OrderBy(p => p.Mz).ToArray();
            var spectrum = new MzSpectrum(peaks.Select(p => p.Mz).ToArray(), peaks.Select(p => p.Height).ToArray(), false);
            return new MsDataScan(spectrum, scan + 1, 2, true, Polarity.Positive, scan, new MzRange(150, 1500), "synthetic", MZAnalyzerType.Orbitrap,
                peaks.Sum(p => p.Height), 10, null, $"scan={scan + 1}", isolationMZ: 500, isolationWidth: 25, dissociationType: DissociationType.HCD);
        }).ToArray();

    private static MslLibraryEntry Entry(params MslFragmentIon[] ions)
    {
        var entry = new MslLibraryEntry { FullSequence = "PEPTIDEK", BaseSequence = "PEPTIDEK", PrecursorMz = 500, ChargeState = 2 };
        entry.MatchedFragmentIons.AddRange(ions);
        return entry;
    }

    /// <summary>The apex is read over three scans in the window holding it, and an ion is read within the fragment tolerance.</summary>
    [Test]
    public void AnAnchorsIonsAreReadAtItsApexInItsWindow()
    {
        var y3 = Ion(ProductType.y, 3, mz: 300);
        var y6 = Ion(ProductType.y, 6, mz: 600);
        var far = Ion(ProductType.y, 7, mz: 700);
        var holding = Window((300, s => s is 9 or 10 or 11 ? 100 : 5000), (600, s => s == 10 ? 150 : 0), (700.1, s => 1000));
        var other = holding.Select(s => new MsDataScan(new MzSpectrum(Array.Empty<double>(), Array.Empty<double>(), false), s.OneBasedScanNumber + 100, 2, true, Polarity.Positive, s.RetentionTime + 0.5,
            new MzRange(150, 1500), "synthetic", MZAnalyzerType.Orbitrap, 0, 10, null, $"scan={s.OneBasedScanNumber + 100}", isolationMZ: 500, isolationWidth: 25, dissociationType: DissociationType.HCD)).ToArray();

        var read = DiaFragmentIntensityModel.ReadAnchor(Entry(y3, y6, far), new RtMinutes(10), [other, holding], 20, true);

        Assert.That(read.Single(r => r.Ion == y3).Observed, Is.EqualTo(300), "three scans of 100");
        Assert.That(read.Single(r => r.Ion == y6).Observed, Is.EqualTo(150));
        Assert.That(read.Single(r => r.Ion == far).Observed, Is.EqualTo(0), "142 ppm off");
    }

    [Test]
    public void AnAnchorWithNoWindowReadsNothing()
    {
        Assert.That(DiaFragmentIntensityModel.ReadAnchor(Entry(Ion(ProductType.y, 3)), new RtMinutes(10), [], 20, true), Is.Empty);
    }

    // Calibration fitting the table: DiaIrtSelfCalibrationTests
}
