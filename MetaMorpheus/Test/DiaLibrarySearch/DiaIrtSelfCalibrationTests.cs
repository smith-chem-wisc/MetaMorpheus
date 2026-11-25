#nullable enable
using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DiaLibrarySearch;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;

namespace Test.DiaLibrarySearch;

/// <summary>
/// A DIA run must be searchable with nothing but its scans and a library. It calibrates its own minutes onto the
/// library's iRT from a wide first pass, and nothing is borrowed from another search engine.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaIrtSelfCalibrationTests
{
    private string _directory = "";

    [OneTimeSetUp]
    public void SetUp()
    {
        _directory = Path.Combine(TestContext.CurrentContext.TestDirectory, "DiaIrtSelfCalibrationTests");
        Directory.CreateDirectory(_directory);
    }

    [OneTimeTearDown]
    public void TearDown()
    {
        Directory.Delete(_directory, true);
    }

    /// <summary>
    /// Calibration scores every fragment whatever the search's fragment rule: with the rule in calibration too, the held-out
    /// run lost 2.1% (two-seed means at a matched paired FDP of 1%). Calibration's parameters are the search's with the rule off,
    /// and nothing else changed; and calibration gives the same result whichever way the search sets the rule.
    /// </summary>
    [Test]
    public void CalibrationIgnoresTheFragmentRule()
    {
        var search = new DiaLibrarySearchParameters(DiaNnFragmentFilter: true, TopFragmentCount: 10);
        var calibrationParameters = DiaIrtSelfCalibration.CalibrationParameters(search);
        Assert.That(calibrationParameters.DiaNnFragmentFilter, Is.False);
        // Candidate spacing, peptidoform evidence and protein support were benchmarked in the main search only: calibration keeps
        // the half width, no rivals and no protein support
        Assert.That(calibrationParameters.CandidateApexSpacingScans, Is.Null);
        Assert.That(calibrationParameters.PeptidoformEvidence, Is.False);
        Assert.That(calibrationParameters.ProteinSupport, Is.False);
        Assert.That(calibrationParameters with { DiaNnFragmentFilter = true, Ms1TolerancePpm = search.Ms1TolerancePpm, DiaNnSignalShare = search.DiaNnSignalShare,
                CandidateApexSpacingScans = search.CandidateApexSpacingScans, PeptidoformEvidence = search.PeptidoformEvidence, ProteinSupport = search.ProteinSupport,
                InterferenceExplainedFragments = search.InterferenceExplainedFragments, ClassifierNetworkLayers = search.ClassifierNetworkLayers,
                ClassifierTreeBlend = search.ClassifierTreeBlend, CandidatePeaksFromSummedCore = search.CandidatePeaksFromSummedCore, NeighbourEvidence = search.NeighbourEvidence, FragmentIntensities = search.FragmentIntensities, IrtComposition = search.IrtComposition,
                MaxToleranceCoElution = search.MaxToleranceCoElution, Ms1TightCorrelation = search.Ms1TightCorrelation, SqrtCoElution = search.SqrtCoElution, Ms1PeakEnvelopeTight = search.Ms1PeakEnvelopeTight, Ms1PeakSignalFraction = search.Ms1PeakSignalFraction, Ms1PeakPoints = search.Ms1PeakPoints, Ms1PeakFeatures = search.Ms1PeakFeatures, Ms1InterpolatedTraces = search.Ms1InterpolatedTraces },
            Is.EqualTo(search), "nothing else changes but calibration's own MS1 tolerance, the signal share, the interference rule, the network, the tree blend, the summed core and the MS1 features it was benchmarked without");

        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var on = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(DiaNnFragmentFilter: true), new CommonParameters());
        var off = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(DiaNnFragmentFilter: false), new CommonParameters());
        Assert.That(on.Anchors, Is.EqualTo(off.Anchors));
    }

    /// <summary>
    /// Targets without a library time (NaN, as mzLib reads a missing one) give calibration no range. The first
    /// pass's window comes from the 2nd and 98th percentiles of the targets' iRTs, and a sort puts NaN below every number,
    /// so a library with a fifth of its targets untimed had a NaN window and no anchors. Calibration reads only the targets
    /// that have a time, so such a library calibrates exactly as one without them.
    /// </summary>
    [Test]
    public void TargetsWithoutALibraryTimeDoNotMoveCalibration()
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var timedLibrary = MslLibrary.Load(run.WriteLibrary(_directory));
        var timed = DiaIrtSelfCalibration.Calibrate(run.Scans, timedLibrary, new DiaLibrarySearchParameters(), new CommonParameters());

        var targets = run.Library.Where(entry => !entry.IsDecoy).ToList();
        run.Library.AddRange(targets.Take(targets.Count / 4).Select((entry, i) => new MslLibraryEntry
        {
            FullSequence = "NANTIME" + entry.FullSequence,
            BaseSequence = "NANTIME" + entry.BaseSequence,
            PrecursorMz = entry.PrecursorMz,
            ChargeState = entry.ChargeState,
            RetentionTime = double.NaN,
            IsDecoy = false,
            MatchedFragmentIons = entry.MatchedFragmentIons,
        }));
        using var untimedLibrary = MslLibrary.Load(run.WriteLibrary(_directory));
        Assert.That(untimedLibrary.TargetCount, Is.EqualTo(timedLibrary.TargetCount + targets.Count / 4));
        var untimed = DiaIrtSelfCalibration.Calibrate(run.Scans, untimedLibrary, new DiaLibrarySearchParameters(), new CommonParameters());

        Assert.That(untimed.IrtHalfWindow, Is.EqualTo(timed.IrtHalfWindow));
        Assert.That(untimed.Anchors, Is.EqualTo(timed.Anchors));
    }

    /// <summary>
    /// Calibration scores without DIA-NN's signal share, whatever the search does: the share was benchmarked in the main search
    /// only (+0.9% / +0.6%), and switching it on in calibration too moved the anchors and cost HF-X 1.5% on seed 0 (58,284
    /// against 59,170). Calibration gives the same anchors whichever way the search sets it.
    /// </summary>
    [Test]
    public void CalibrationIgnoresTheSignalShare()
    {
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters(DiaNnSignalShare: true)).DiaNnSignalShare, Is.False);
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters()).MaxToleranceCoElution, Is.False, "MaxToleranceCoElution: calibration was benchmarked without it");
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters()).Ms1TightCorrelation, Is.False, "Ms1TightCorrelation: calibration was benchmarked without it");
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters()).SqrtCoElution, Is.False, "SqrtCoElution: calibration was benchmarked without it");
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters()).Ms1PeakEnvelopeTight, Is.False, "Ms1PeakEnvelopeTight: calibration was benchmarked without it");
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters()).Ms1PeakSignalFraction, Is.False, "Ms1PeakSignalFraction: calibration was benchmarked without it");
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters()).Ms1PeakPoints, Is.False, "Ms1PeakPoints: calibration was benchmarked without it");
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters()).Ms1PeakFeatures, Is.False, "Ms1PeakFeatures: calibration was benchmarked without it");
        var ms1Options = DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters(Ms1InterpolatedTraces: true, Ms1CorrelationSummedCore: true));
        Assert.That(ms1Options.Ms1InterpolatedTraces, Is.False, "Ms1InterpolatedTraces: calibration was benchmarked without it");
        Assert.That(ms1Options.Ms1CorrelationSummedCore, Is.False, "Ms1CorrelationSummedCore: calibration was benchmarked without it");

        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var on = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(DiaNnSignalShare: true), new CommonParameters());
        var off = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(DiaNnSignalShare: false), new CommonParameters());
        Assert.That(on.Anchors, Is.EqualTo(off.Anchors));
    }

    /// <summary>
    /// Calibration removes interference at 4 explained fragments, whatever the search uses: the main search's 3 was
    /// benchmarked with calibration at 4. Calibration gives the same anchors whichever rule the search sets.
    /// </summary>
    [Test]
    public void CalibrationKeepsItsOwnInterferenceRule()
    {
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters(InterferenceExplainedFragments: 3)).InterferenceExplainedFragments,
            Is.EqualTo(DiaIrtSelfCalibration.CalibrationInterferenceExplainedFragments));
        Assert.That(DiaIrtSelfCalibration.CalibrationInterferenceExplainedFragments, Is.EqualTo(4));
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters(InterferenceSameMzOnly: true)).InterferenceSameMzOnly, Is.True,
            "the same-m/z pairing was benchmarked with calibration pairing the same way");

        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, noisePeaksPerScan: 3000);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var three = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(InterferenceExplainedFragments: 3), new CommonParameters());
        var six = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(InterferenceExplainedFragments: 6), new CommonParameters());
        Assert.That(three.Anchors, Is.EqualTo(six.Anchors));
    }

    /// <summary>
    /// Calibration's network keeps the rescorer's 25-20-15-10-5 (DIA-NN 2020's), whatever the search uses: the main search's
    /// 64-32-16 was benchmarked with calibration at 25-20-15-10-5. Calibration gives the same anchors whichever network the
    /// search sets.
    /// </summary>
    [Test]
    public void CalibrationKeepsItsOwnNetwork()
    {
        Assert.That(DiaIrtSelfCalibration.CalibrationNetworkLayers, Is.EqualTo(new[] { 25, 20, 15, 10, 5 }));
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters()).EffectiveClassifierNetworkLayers,
            Is.EqualTo(DiaIrtSelfCalibration.CalibrationNetworkLayers));
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters(ClassifierNetworkLayers: [16, 8])).EffectiveClassifierNetworkLayers,
            Is.EqualTo(DiaIrtSelfCalibration.CalibrationNetworkLayers));

        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var standard = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters());
        var narrow = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(ClassifierNetworkLayers: [16, 8]), new CommonParameters());
        Assert.That(narrow.Anchors, Is.EqualTo(standard.Anchors));
    }

    /// <summary>
    /// Calibration also fits the run's MS1 offset from its anchors' raw MS1 errors; the synthetic run has none, so the
    /// fitted offset is near 0. Calibration's own passes search without an offset.
    /// </summary>
    [Test]
    public void CalibrationFitsTheRunsMs1Offset()
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, withMs1: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters());

        Assert.That(calibration.Ms1Offset, Is.Not.Null);
        Assert.That(Math.Abs(calibration.Ms1Offset!.OffsetPpm(5)), Is.LessThan(1));
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters(Ms1Offset: Ms1OffsetModel.Constant(3))).Ms1Offset, Is.Null);
        // Calibration reads MS1 at its own fixed tolerance (the offset is not known yet), whatever the search's
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters(Ms1TolerancePpm: 5)).Ms1TolerancePpm,
            Is.EqualTo(DiaIrtSelfCalibration.CalibrationMs1TolerancePpm));
        // ...unless the search is set wider: a user who knows the instrument is off widens calibration's read too
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters(Ms1TolerancePpm: 30)).Ms1TolerancePpm, Is.EqualTo(30));
        Assert.That(calibration.Ms1OffsetAnchors, Is.GreaterThan(calibration.AnchorCount / 2), "most anchors have an M0 the offset is fitted from");
        // Applying a calibration sets the search's iRT window, MS1 offset, intensity table and iRT composition fit together
        var applied = calibration.ApplyTo(new DiaLibrarySearchParameters());
        Assert.That(applied.IrtHalfWindow, Is.EqualTo(calibration.IrtHalfWindow));
        Assert.That(applied.Ms1Offset, Is.SameAs(calibration.Ms1Offset));
        Assert.That(calibration.FragmentIntensities, Is.Not.Null);
        Assert.That(applied.FragmentIntensities, Is.SameAs(calibration.FragmentIntensities));
        Assert.That(calibration.IrtComposition, Is.Not.Null);
        Assert.That(applied.IrtComposition, Is.SameAs(calibration.IrtComposition));
    }

    /// <summary>
    /// A run whose MS1 error is past calibration's 10 ppm read: almost no anchor found its M0 there, so the offset was fitted
    /// from a couple of chance peaks, near 0 (12 ppm fitted 0.07), and nothing said so. The offset is then read from each
    /// anchor's nearest MS1 scan at the 20 ppm mass calibration reads with, and it comes out right, resting on most anchors.
    /// A search set to 30 ppm reads that far in calibration too, and its own read is enough.
    /// </summary>
    [TestCase(12, 5)]
    [TestCase(-15, 5)]
    [TestCase(12, 30)]
    public void AnMs1OffsetPastCalibrationsReadIsStillFitted(double ppm, double searchTolerancePpm)
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, withMs1: true, ppmOffset: ppm);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(Ms1TolerancePpm: searchTolerancePpm), new CommonParameters());

        Assert.That(calibration.Ms1Offset, Is.Not.Null);
        Assert.That(calibration.Ms1Offset!.OffsetPpm(5), Is.EqualTo(ppm).Within(1));
        Assert.That(calibration.Ms1OffsetAnchors, Is.GreaterThan(calibration.AnchorCount / 2));
    }

    /// <summary>The library halves one class's prediction after the run was made, so the run reads that class at twice what the library says, and calibration's table says so.</summary>
    [Test]
    public void CalibrationFitsTheIntensityTableFromItsAnchors()
    {
        var run = SyntheticDiaRun.Build(1000, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 60);
        foreach (var ion in run.Library.SelectMany(e => e.MatchedFragmentIons).Where(f => f.FragmentNumber == 5))
            ion.Intensity /= 2;
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters());

        Assert.That(calibration.FragmentIntensities, Is.Not.Null);
        double ratio = DiaFragmentIntensityModelTests.Ratio(calibration.FragmentIntensities!,
            DiaFragmentIntensityModelTests.Ion(ProductType.y, 5), DiaFragmentIntensityModelTests.Ion(ProductType.y, 6));
        Assert.That(ratio, Is.EqualTo(2).Within(0.3));

        var oneRound = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(), rounds: 1);
        Assert.That(oneRound.FragmentIntensities, Is.Null, "as with the MS1 offset, only a later round fits one");
    }

    /// <summary>The library is moved 5 iRT earlier per W after the run was made, so the run reads those precursors 5 iRT later than it says, and calibration's fit says so.</summary>
    [Test]
    public void CalibrationFitsTheIrtCompositionFromItsAnchors()
    {
        var run = SyntheticDiaRun.Build(1000, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 60);
        foreach (var entry in run.Library.Where(e => e.BaseSequence.Contains('W')))
            entry.RetentionTime -= 5 * entry.BaseSequence.Count(c => c == 'W');
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters());

        Assert.That(calibration.IrtComposition, Is.Not.Null);
        Assert.That(calibration.IrtComposition!.Offset("AAAAWAAAAK") - calibration.IrtComposition.Offset("AAAAAAAAAK"), Is.EqualTo(5).Within(1.5));

        var oneRound = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(), rounds: 1);
        Assert.That(oneRound.IrtComposition, Is.Null, "as with the MS1 offset, only a later round fits one");
    }

    /// <summary>
    /// Calibration searches with the settings it was benchmarked with: an intensity table or iRT composition fit passed in is
    /// dropped (calibration fits its own), the main search's tree blend, summed-core candidates and neighbour evidence are
    /// off, and it trains as many networks as the main search, as the A/B ran.
    /// </summary>
    [Test]
    public void CalibrationLeavesOutWhatOnlyTheMainSearchUses()
    {
        var calibrationOf = DiaIrtSelfCalibration.CalibrationParameters;
        Assert.That(calibrationOf(new DiaLibrarySearchParameters(FragmentIntensities: DiaFragmentIntensityModel.Fit([]))).FragmentIntensities, Is.Null);
        Assert.That(calibrationOf(new DiaLibrarySearchParameters(IrtComposition: DiaIrtCompositionModel.Fit([]))).IrtComposition, Is.Null);
        Assert.That(calibrationOf(new DiaLibrarySearchParameters(ClassifierTreeBlend: true)).ClassifierTreeBlend, Is.False);
        Assert.That(calibrationOf(new DiaLibrarySearchParameters(CandidatePeaksFromSummedCore: true)).CandidatePeaksFromSummedCore, Is.False);
        Assert.That(calibrationOf(new DiaLibrarySearchParameters(NeighbourEvidence: true)).NeighbourEvidence, Is.False);
        Assert.That(calibrationOf(new DiaLibrarySearchParameters()).ClassifierNetworkMembers, Is.EqualTo(24));
    }

    /// <summary>The first pass's confident identifications recover the run's true, nonlinear RT(iRT).</summary>
    [Test]
    public void TheFirstPassRecoversTheRunsTrueCurve()
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters());

        for (double irt = -10; irt <= 110; irt += 10)
        {
            var observed = calibration.Model.ToIrt(new RtMinutes(SyntheticDiaRun.TrueRtMinutes(irt)));
            Assert.That(observed.Value, Is.EqualTo(irt).Within(2.0), $"at library iRT {irt}");
        }
        Assert.That(calibration.IrtHalfWindow, Is.LessThan(20), "a calibrated run needs a narrower window than the first pass");
        Assert.That(calibration.AnchorCount, Is.GreaterThanOrEqualTo(50));
        // The anchors are kept, so a bad fit can be traced to the identifications it was fitted on
        Assert.That(calibration.Anchors, Has.Count.EqualTo(calibration.AnchorCount));
        Assert.That(calibration.Anchors.Count(a => Math.Abs(calibration.Model.ToIrt(a.ApexRt).Value - a.LibraryIrt.Value) < 3), Is.GreaterThan(0.9 * calibration.AnchorCount));
    }

    /// <summary>Two passes find what one pass with a borrowed map found, with no map supplied at all.</summary>
    [Test]
    public void ASelfCalibratedSearchFindsThePlantedPrecursors()
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters());
        var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, calibration.Model,
            new DiaLibrarySearchParameters(IrtHalfWindow: calibration.IrtHalfWindow), new CommonParameters(), [], []).Run();

        int found = results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence));
        Assert.That(found, Is.GreaterThanOrEqualTo((int)Math.Ceiling(0.9 * run.PlantedSequences.Count)));
    }

    /// <summary>
    /// The first pass needs only enough targets to yield its anchors. A large library is sampled down to about
    /// <see cref="DiaIrtSelfCalibration.FirstPassTargetCount"/> targets; a small one is searched whole.
    /// </summary>
    [TestCase(300, 1)]
    [TestCase(10_000, 1)]
    [TestCase(10_001, 2)]
    [TestCase(38_987, 4)]
    public void TheFirstPassSamplesALargeLibrary(int targets, int stride)
    {
        Assert.That(DiaIrtSelfCalibration.FirstPassStride(targets), Is.EqualTo(stride));
    }

    /// <summary>
    /// A whole-proteome library is sparse: few of its precursors are in any run. A sample then holds too few for any to
    /// reach 1% q, since (D+1)/T needs a hundred targets before the first decoy. The first pass samples twice as many
    /// precursors until it finds enough anchors, up to the whole library.
    /// </summary>
    [Test]
    public void ASparseLibraryIsSampledMoreUntilTheFirstPassFindsAnchors()
    {
        // 3000 targets of which about 5% elute: a 500-target sample holds about 25
        var run = SyntheticDiaRun.Build(3000, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 20) == 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(),
            firstPassTargetCount: 500);

        Assert.That(calibration.AnchorCount, Is.GreaterThanOrEqualTo(100));
        for (double irt = -10; irt <= 110; irt += 20)
            Assert.That(calibration.Model.ToIrt(new RtMinutes(SyntheticDiaRun.TrueRtMinutes(irt))).Value, Is.EqualTo(irt).Within(2.0), $"at library iRT {irt}");
    }

    /// <summary>
    /// A run can spend minutes before anything elutes, and a short gradient elutes only part of a whole-proteome library's
    /// iRT range. PXD071205's 5-min direct-injection Astral runs do both: peptides elute from about 2.5 of 6.8 min. The first
    /// pass lays the library's iRT range over the whole MS2 time span, so its provisional line then puts every eluting
    /// precursor outside its window, and it finds nothing. When the first pass finds too few anchors to fit, even over the
    /// whole library, calibration searches once more with a window over the library's whole iRT range.
    /// </summary>
    [Test]
    public void ElutionAfterDeadTimeIsStillCalibrated()
    {
        // Peptides below iRT 30 elute after 10 min of dead time in a 20-min run. The provisional line places them 40 to 68
        // iRT from where they elute, beyond the first pass's half-window of about 31. 600 targets leave about 160 eluting:
        // (D+1)/T reaches 1% q only past 100 targets.
        const double deadMinutes = 10;
        var run = SyntheticDiaRun.Build(600, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0 && entry.RetentionTime < 30,
            deadMinutes: deadMinutes);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters());

        for (double irt = -10; irt <= 20; irt += 10)
            Assert.That(calibration.Model.ToIrt(new RtMinutes(deadMinutes + SyntheticDiaRun.TrueRtMinutes(irt))).Value, Is.EqualTo(irt).Within(2.0),
                $"at library iRT {irt}");
    }

    /// <summary>
    /// A sample where the provisional window finds no anchor at all is searched over the library's whole iRT range too, and
    /// calibration continues there when that finds enough. A run with dead time then no longer searches the whole library
    /// with a window that cannot work before it widens. On PXD071205's 5-min runs that cost about 38 of calibration's 47 min.
    /// Here 3000 targets sampled 500 at a time (stride 6): widening only after the whole library takes 6 or more searches,
    /// widening at the first empty sample takes 4.
    /// </summary>
    [Test]
    public void AnEmptyProvisionalWindowWidensAtOnce()
    {
        const double deadMinutes = 10;
        var run = SyntheticDiaRun.Build(3000, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0 && entry.RetentionTime < 30,
            deadMinutes: deadMinutes);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(),
            firstPassTargetCount: 500);

        Assert.That(calibration.FirstPassSearches, Is.LessThanOrEqualTo(4));
        for (double irt = -10; irt <= 20; irt += 10)
            Assert.That(calibration.Model.ToIrt(new RtMinutes(deadMinutes + SyntheticDiaRun.TrueRtMinutes(irt))).Value, Is.EqualTo(irt).Within(2.0),
                $"at library iRT {irt}");
    }

    /// <summary>
    /// A whole-range search can find anchors where the provisional window found none for a reason other than a wrong line: it
    /// searches more of a small sample's precursors, enough to reach 1% q (on the HeLa runs it did, and switching cost
    /// held-out half its anchors). So the switch asks the anchors: the line is misplaced only when most of them sit outside
    /// the window the first pass used.
    /// </summary>
    [TestCase(0.0, false)]
    [TestCase(0.4, false)]
    [TestCase(0.5, false)]
    [TestCase(0.6, true)]
    [TestCase(1.0, true)]
    public void WideAnchorsSayTheProvisionalLineIsMisplacedWhenMostLieOutsideItsWindow(double outsideFraction, bool misplaced)
    {
        var line = Chromatography.RetentionTimeCalibration.IrtCalibration.Line((new RtMinutes(0), new Irt(0)), (new RtMinutes(10), new Irt(100)));
        const double halfWindow = 10;
        int outside = (int)Math.Round(outsideFraction * 50);
        // Each anchor's apex maps to iRT 10 * minutes on the line; library iRT 5 away is inside the window, 40 away is outside
        var anchors = Enumerable.Range(0, 50)
            .Select(i => (new RtMinutes(i * 0.2), new Irt(i * 2.0 + (i < outside ? 40 : 5))))
            .ToList();

        Assert.That(DiaIrtSelfCalibration.LineIsMisplaced(anchors, line, halfWindow), Is.EqualTo(misplaced));
    }

    /// <summary>A run whose provisional line is right takes its first sample, and is never searched over the whole range.</summary>
    [Test]
    public void ARunThatCalibratesOnItsFirstSampleSearchesOnce()
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(), desiredAnchors: 20);

        Assert.That(calibration.FirstPassSearches, Is.EqualTo(1));
    }

    /// <summary>Calibration can be repeated on a disjoint first-pass sample; each sample calibrates the run on its own.</summary>
    [TestCase(1)]
    [TestCase(3)]
    public void AnotherFirstPassSampleCalibratesToo(int offset)
    {
        var run = SyntheticDiaRun.Build(3000, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) == 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(),
            firstPassTargetCount: 500, sampleOffset: offset);

        Assert.That(calibration.AnchorCount, Is.GreaterThanOrEqualTo(100));
        for (double irt = -10; irt <= 110; irt += 20)
            Assert.That(calibration.Model.ToIrt(new RtMinutes(SyntheticDiaRun.TrueRtMinutes(irt))).Value, Is.EqualTo(irt).Within(2.0), $"at library iRT {irt}");
        Assert.Throws<ArgumentOutOfRangeException>(() => DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(),
            new CommonParameters(), sampleOffset: -1));
    }

    /// <summary>
    /// The first pass widens until it has the anchors asked for. A run whose first pass stopped at 843 anchors needed a second
    /// search to fix its RT model (PXD022589, +18%); asking for more anchors up front is the cheaper fix to try.
    /// </summary>
    [Test]
    public void AskingForMoreAnchorsWidensTheFirstPassFurther()
    {
        var run = SyntheticDiaRun.Build(3000, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) == 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var few = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(),
            firstPassTargetCount: 500, desiredAnchors: 100);
        var many = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(),
            firstPassTargetCount: 500, desiredAnchors: 600);

        Assert.That(few.AnchorCount, Is.GreaterThanOrEqualTo(100).And.LessThan(600));
        Assert.That(many.AnchorCount, Is.GreaterThanOrEqualTo(600));
    }

    /// <summary>
    /// The search window is a multiple of the calibration's residual SD (at least 5 iRT). On PXD005573, every change that
    /// narrowed it cost identifications, so the multiple is a parameter to measure rather than a constant.
    /// </summary>
    [TestCase(4.0)]
    [TestCase(5.0)]
    public void TheWindowIsTheChosenMultipleOfTheResidualSd(double sds)
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(), windowSds: sds);

        Assert.That(calibration.IrtHalfWindow, Is.EqualTo(Math.Max(5, sds * calibration.Model.ResidualSd)).Within(1e-9));
    }

    /// <summary>
    /// Calibration passes pick anchors with the linear discriminant, whatever model the main search uses: anchors only need to
    /// be confidently right. On the whole-proteome library this cut calibration from 12:18 to 3:51 (PXD005573 1 h) and gave
    /// slightly more precursors at matched FDP.
    /// </summary>
    [Test]
    public void CalibrationUsesTheLinearModelWhateverTheSearchUses()
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var network = DiaIrtSelfCalibration.Calibrate(run.Scans, library,
            new DiaLibrarySearchParameters(ClassifierModel: StatisticalModels.RescoreModel.NeuralNetworkEnsemble), new CommonParameters());
        var linear = DiaIrtSelfCalibration.Calibrate(run.Scans, library,
            new DiaLibrarySearchParameters(ClassifierModel: StatisticalModels.RescoreModel.LinearDiscriminant), new CommonParameters());

        Assert.That(network.AnchorCount, Is.EqualTo(linear.AnchorCount));
        Assert.That(network.IrtHalfWindow, Is.EqualTo(linear.IrtHalfWindow));
        Assert.That(network.Model.ResidualSd, Is.EqualTo(linear.Model.ResidualSd));
    }

    /// <summary>
    /// Iterative calibration (DIA-NN calibrates in rounds): each later round searches the same sample again with the previous
    /// round's fitted curve and window instead of the provisional straight line, so anchors come from the whole gradient
    /// rather than where the straight line happened to be right. The curve stays accurate, and the window is the last round's.
    /// </summary>
    [Test]
    public void ALaterCalibrationRoundSearchesWithTheFittedCurve()
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var one = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(), rounds: 1);
        var two = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(), rounds: 2);

        Assert.That(two.Rounds, Is.EqualTo(2));
        Assert.That(one.Rounds, Is.EqualTo(1));
        Assert.That(two.AnchorCount, Is.GreaterThanOrEqualTo(one.AnchorCount * 0.9));
        for (double irt = -10; irt <= 110; irt += 10)
            Assert.That(two.Model.ToIrt(new RtMinutes(SyntheticDiaRun.TrueRtMinutes(irt))).Value, Is.EqualTo(irt).Within(2.0), $"at library iRT {irt}");
    }

    /// <summary>
    /// The second pass (DIA-NN refits after its first search): the main search's confident targets refit the run's RT->iRT map.
    /// Only targets at 1% count. A far-off target above 1%, or any decoy, must not pull the fit.
    /// </summary>
    [Test]
    public void ARefitFromTheMainSearchUsesOnlyItsConfidentTargets()
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters());
        var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, calibration.Model,
            new DiaLibrarySearchParameters(IrtHalfWindow: calibration.IrtHalfWindow), new CommonParameters(), [], []).Run();
        int confident = results.Matches.Count(m => !m.IsDecoy && m.QValue <= DiaIrtSelfCalibration.AnchorQValue);
        // Poison: wrong apexes far from the curve, as a weak target and as a decoy
        var poisoned = results.Matches.Concat(Enumerable.Range(0, 50).SelectMany(i => new[]
        {
            new DiaPrecursorMatch(100_000 + i, "POISONK", 2, 500, false, new Irt(100), new RtMinutes(1), new Irt(0), 1, [1.0], QValue: 0.5),
            new DiaPrecursorMatch(200_000 + i, "KNOSIOP", 2, 500, true, new Irt(100), new RtMinutes(1), new Irt(0), 99, [1.0]),
        })).ToList();

        var refined = DiaIrtSelfCalibration.Refine(calibration, poisoned, library);

        Assert.That(refined.AnchorCount, Is.EqualTo(Math.Min(confident, DiaIrtSelfCalibration.MaximumRefinementAnchors)));
        for (double irt = -10; irt <= 110; irt += 10)
            Assert.That(refined.Model.ToIrt(new RtMinutes(SyntheticDiaRun.TrueRtMinutes(irt))).Value, Is.EqualTo(irt).Within(2.0), $"at library iRT {irt}");
    }

    /// <summary>
    /// A refit replaces the curve and its window, not the rest of the calibration: applied on top of the first calibration,
    /// it keeps that calibration's MS1 offset and intensity table, and its iRT composition fit is refitted against the new
    /// curve. A refit with too few targets says so, not that the first pass found too few.
    /// </summary>
    [Test]
    public void ARefitKeepsTheCalibrationsOtherModels()
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, withMs1: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var calibration = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters());
        Assert.That(calibration.Ms1Offset, Is.Not.Null);
        Assert.That(calibration.FragmentIntensities, Is.Not.Null);
        Assert.That(calibration.IrtComposition, Is.Not.Null);
        var first = calibration.ApplyTo(new DiaLibrarySearchParameters());
        var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, calibration.Model, first, new CommonParameters(), [], []).Run();

        var refined = DiaIrtSelfCalibration.Refine(calibration, results.Matches, library);
        var second = refined.ApplyTo(first);

        Assert.That(second.IrtHalfWindow, Is.EqualTo(refined.IrtHalfWindow));
        Assert.That(second.Ms1Offset, Is.SameAs(calibration.Ms1Offset));
        Assert.That(second.FragmentIntensities, Is.SameAs(calibration.FragmentIntensities));
        Assert.That(second.IrtComposition, Is.Not.Null.And.Not.SameAs(calibration.IrtComposition), "refitted against the new curve");
        Assert.That(refined.Ms1OffsetAnchors, Is.EqualTo(calibration.Ms1OffsetAnchors));
        var oneRound = DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters(), rounds: 1);
        Assert.That(DiaIrtSelfCalibration.Refine(oneRound, results.Matches, library).IrtComposition, Is.Null, "a calibration without a composition fit gets none");
        var e = Assert.Throws<MetaMorpheusException>(() => DiaIrtSelfCalibration.Refine(calibration, [], library));
        Assert.That(e!.Message, Does.Contain("main search").And.Not.Contain("first pass"));
    }

    /// <summary>With nothing to find there is nothing to calibrate on, and that is an error, not a guess.</summary>
    [Test]
    public void ARunWithNothingToFindCannotBeCalibrated()
    {
        var run = SyntheticDiaRun.Build(200, _ => false);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var e = Assert.Throws<MetaMorpheusException>(() =>
            DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(), new CommonParameters()));
        Assert.That(e!.Message, Does.Contain("calibrat"));
    }

    /// <summary>
    /// What cannot be calibrated fails before any search, and says why: no calibration rounds, a run with no DIA MS2 scans, and
    /// a library whose targets have no library time to set the first pass's range.
    /// </summary>
    [Test]
    public void WhatCannotBeCalibratedFailsBeforeAnySearch()
    {
        var run = SyntheticDiaRun.Build(50, _ => false, withMs1: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        Assert.Throws<ArgumentOutOfRangeException>(() => DiaIrtSelfCalibration.Calibrate(run.Scans, library, new DiaLibrarySearchParameters(),
            new CommonParameters(), rounds: 0));

        var ms1Only = run.Scans.Where(s => s.MsnOrder == 1).ToArray();
        Assert.That(ms1Only, Is.Not.Empty);
        var noMs2 = Assert.Throws<MetaMorpheusException>(() =>
            DiaIrtSelfCalibration.Calibrate(ms1Only, library, new DiaLibrarySearchParameters(), new CommonParameters()));
        Assert.That(noMs2!.Message, Does.Contain("no DIA MS2 scans"));

        var untimedTargets = run.Library.Where(entry => !entry.IsDecoy).Select(entry => new MslLibraryEntry
        {
            FullSequence = entry.FullSequence,
            BaseSequence = entry.BaseSequence,
            PrecursorMz = entry.PrecursorMz,
            ChargeState = entry.ChargeState,
            RetentionTime = double.NaN,
            IsDecoy = false,
            MatchedFragmentIons = entry.MatchedFragmentIons,
        }).ToList();
        run.Library.RemoveAll(entry => !entry.IsDecoy);
        run.Library.AddRange(untimedTargets);
        using var untimedLibrary = MslLibrary.Load(run.WriteLibrary(_directory));
        Assert.That(untimedLibrary.TargetCount, Is.GreaterThan(0));
        var untimed = Assert.Throws<MetaMorpheusException>(() =>
            DiaIrtSelfCalibration.Calibrate(run.Scans, untimedLibrary, new DiaLibrarySearchParameters(), new CommonParameters()));
        Assert.That(untimed!.Message, Does.Contain("no target precursors with a library time"));
    }
}
