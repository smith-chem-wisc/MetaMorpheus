#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Security.Cryptography;
using EngineLayer;
using Chromatography.RetentionTimeCalibration;
using EngineLayer.DiaLibrarySearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;

namespace Test.DiaLibrarySearch;

/// <summary>
/// The DIA library search's permanent contract, run on <see cref="SyntheticDiaRun"/>, where the answer is known. Every
/// later milestone deepens the engine; none may break these.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaLibrarySearchEngineTests
{
    private string _directory = "";

    /// <summary>
    /// A linear run-RT to iRT map through the true curve's endpoints (iRT -20 at 0.5 min, iRT 120 at 9.936 min). The
    /// true curve bows away from it by up to about 12 iRT, inside the engine's default iRT window.
    /// </summary>
    private static readonly IrtCalibrationModel EndpointMap = IrtCalibration.Line(
        (new RtMinutes(SyntheticDiaRun.TrueRtMinutes(-20)), new Irt(-20)),
        (new RtMinutes(SyntheticDiaRun.TrueRtMinutes(120)), new Irt(120)));

    [OneTimeSetUp]
    public void SetUp()
    {
        _directory = Path.Combine(TestContext.CurrentContext.TestDirectory, "DiaLibrarySearchEngineTests");
        Directory.CreateDirectory(_directory);
    }

    [OneTimeTearDown]
    public void TearDown()
    {
        Directory.Delete(_directory, true);
    }

    private DiaLibrarySearchResults Search(SyntheticDiaRun run, IrtCalibrationModel map, out string libraryPath)
    {
        libraryPath = run.WriteLibrary(_directory);
        using var library = MslLibrary.Load(libraryPath);
        var engine = new DiaLibrarySearchEngine(run.Scans, library, map, new DiaLibrarySearchParameters(),
            new CommonParameters(), [], []);
        return (DiaLibrarySearchResults)engine.Run();
    }

    private DiaLibrarySearchResults Search(SyntheticDiaRun run) => Search(run, EndpointMap, out _);

    private static List<DiaPrecursorMatch> Search(MslLibrary library, IEnumerable<MsDataScan> scans) =>
        ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(scans.ToArray(), library, EndpointMap, new DiaLibrarySearchParameters(),
            new CommonParameters(), [], []).Run()).Matches.OrderBy(m => m.PrecursorIndex).ToList();

    /// <summary>Three of every four targets elute. At least 90% of those are reported at 1% FDR.</summary>
    [Test]
    public void PlantedPrecursorsAreFound()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);

        var results = Search(run);

        var found = results.Matches.Where(m => !m.IsDecoy && m.QValue <= 0.01).Select(m => m.FullSequence).ToHashSet();
        int planted = run.PlantedSequences.Count;
        Assert.That(planted, Is.GreaterThan(100), "the fixture must plant enough precursors to mean something");
        Assert.That(found.Count(run.PlantedSequences.Contains), Is.GreaterThanOrEqualTo((int)Math.Ceiling(0.9 * planted)));
        Assert.That(found.Count(sequence => !run.PlantedSequences.Contains(sequence)), Is.LessThanOrEqualTo(Math.Max(1, found.Count / 100)),
            "a target that never eluted is a false discovery; at 1% FDR there can be at most about 1% of them");
    }

    /// <summary>A run of noise alone has nothing to find.</summary>
    [Test]
    public void NullRunReportsNothing()
    {
        var run = SyntheticDiaRun.Build(200, _ => false);

        var results = Search(run);

        Assert.That(results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01), Is.LessThanOrEqualTo(Math.Max(1, results.Matches.Count / 100)));
    }

    /// <summary>When only the decoys elute, no target may pass: every real signal belongs to a decoy.</summary>
    [Test]
    public void ARunOfDecoysReportsNoTargets()
    {
        var run = SyntheticDiaRun.Build(200, entry => entry.IsDecoy);

        var results = Search(run);

        Assert.That(results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01), Is.Zero);
    }

    /// <summary>A library with no decoys gives the search nothing to estimate an FDR from, so it is refused, not searched.</summary>
    [Test]
    public void ALibraryWithoutDecoysThrows()
    {
        var run = SyntheticDiaRun.Build(50, entry => true, withDecoys: false);

        var e = Assert.Throws<MetaMorpheusException>(() => Search(run));
        Assert.That(e!.Message, Does.Contain("decoy"));
    }

    /// <summary>
    /// A library's decoys can all go unscored: here none elutes and the run has no noise for one to match. D is then 0,
    /// so (D+1)/T, made monotone, gives every target 1/T. Five targets cannot pass 1% FDR however well they score.
    /// </summary>
    [TestCase(5)]
    [TestCase(12)]
    public void WithNoDecoyScoredEveryTargetGetsOneOverT(int targetCount)
    {
        var run = SyntheticDiaRun.Build(targetCount, entry => !entry.IsDecoy, noisePeaksPerScan: 0);

        var results = Search(run);

        Assert.That(run.Library.Count(entry => entry.IsDecoy), Is.EqualTo(targetCount), "the library holds decoys");
        Assert.That(results.DecoyCount, Is.Zero, "but the search scored none");
        Assert.That(results.TargetCount, Is.EqualTo(targetCount));
        Assert.That(results.Matches.Select(m => m.QValue), Is.All.EqualTo(1.0 / targetCount));
    }

    /// <summary>
    /// Decoys can also go unscored by lying outside every searched window, in m/z or in iRT, so that the search never
    /// reaches them. Left where they were, the same decoys score on this run's noise (dense enough, at 3000 peaks a scan,
    /// for a decoy's candidate peak to form).
    /// </summary>
    [TestCase(true)]
    [TestCase(false)]
    public void DecoysOutsideEverySearchedWindowLeaveEveryTargetOneOverT(bool outsideInMz)
    {
        var run = SyntheticDiaRun.Build(20, entry => !entry.IsDecoy, noisePeaksPerScan: 3000);
        Assert.That(Search(run).DecoyCount, Is.Positive, "left where they were, the decoys score");
        foreach (var decoy in run.Library.Where(entry => entry.IsDecoy))
        {
            if (outsideInMz)
                decoy.PrecursorMz = SyntheticDiaRun.FirstWindowLowMz - 100;
            else
                decoy.RetentionTime = 10000;
        }

        var results = Search(run);

        Assert.That(results.DecoyCount, Is.Zero);
        Assert.That(results.TargetCount, Is.EqualTo(20));
        Assert.That(results.Matches.Select(m => m.QValue), Is.All.EqualTo(1.0 / 20));
    }

    /// <summary>
    /// Enough targets can pass on (0+1)/T alone, so a search that scored targets but no decoy says so. An ordinary search
    /// does not, and neither does one that scored nothing, whose results claim nothing.
    /// </summary>
    [Test]
    [NonParallelizable]
    public void OnlyASearchThatScoredTargetsButNoDecoyWarns()
    {
        var warnings = new List<string>();
        EventHandler<StringEventArgs> onWarn = (_, e) => warnings.Add(e.S);
        var withoutDecoyScores = SyntheticDiaRun.Build(5, entry => !entry.IsDecoy, noisePeaksPerScan: 0);
        var ordinary = SyntheticDiaRun.Build(50, entry => !entry.IsDecoy, noisePeaksPerScan: 3000);
        using var withoutDecoyScoresLibrary = MslLibrary.Load(withoutDecoyScores.WriteLibrary(_directory));
        using var ordinaryLibrary = MslLibrary.Load(ordinary.WriteLibrary(_directory));
        MetaMorpheusEngine.WarnHandler += onWarn;
        try
        {
            Search(ordinaryLibrary, ordinary.Scans);
            Search(ordinaryLibrary, []);
            Assert.That(warnings.Where(w => w.Contains("No decoy precursor was scored")), Is.Empty);

            Search(withoutDecoyScoresLibrary, withoutDecoyScores.Scans);
            Assert.That(warnings.Count(w => w.Contains("No decoy precursor was scored")), Is.EqualTo(1));
        }
        finally
        {
            MetaMorpheusEngine.WarnHandler -= onWarn;
        }
    }

    /// <summary>
    /// A sampled search, as calibration's first pass uses, scores only precursors whose library index is a multiple of the
    /// stride. It samples targets and decoys alike, so the target-decoy competition stays fair.
    /// </summary>
    [Test]
    public void ASampledSearchScoresOnlyEveryStrideThPrecursor()
    {
        // Decoys elute too: an entry with no signal at all is never scored, and the test needs both classes
        var run = SyntheticDiaRun.Build(200, entry => SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var full = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(), new CommonParameters(), [], []).Run();
        var sampled = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(PrecursorSampleStride: 3), new CommonParameters(), [], []).Run();

        Assert.That(sampled.Matches.All(m => m.PrecursorIndex % 3 == 0));
        Assert.That(sampled.TargetCount, Is.EqualTo(full.Matches.Count(m => !m.IsDecoy && m.PrecursorIndex % 3 == 0)));
        Assert.That(sampled.DecoyCount, Is.EqualTo(full.Matches.Count(m => m.IsDecoy && m.PrecursorIndex % 3 == 0)));
        Assert.That(sampled.DecoyCount, Is.GreaterThan(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DiaLibrarySearchParameters(PrecursorSampleStride: 0));
    }

    /// <summary>
    /// An offset picks which of the stride's samples is scored, so calibration can be repeated on disjoint samples: a
    /// different first-pass sample cost HF-X 20% of its precursors, and the spread over samples is what to measure.
    /// </summary>
    [Test]
    public void AnOffsetPicksADisjointSample()
    {
        var run = SyntheticDiaRun.Build(200, entry => SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var samples = Enumerable.Range(0, 3).Select(offset => ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(PrecursorSampleStride: 3) { PrecursorSampleOffset = offset }, new CommonParameters(), [], []).Run()).Matches).ToList();

        for (int offset = 0; offset < 3; offset++)
            Assert.That(samples[offset].All(m => m.PrecursorIndex % 3 == offset), $"offset {offset}");
        Assert.That(samples.Sum(s => s.Count), Is.EqualTo(((DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(), new CommonParameters(), [], []).Run()).Matches.Count));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DiaLibrarySearchParameters(PrecursorSampleStride: 3) { PrecursorSampleOffset = 3 });
    }

    /// <summary>
    /// The classifier's network can train twice, as DIA-NN's does: the second pass learns from the candidate peaks the first
    /// network picked. The setting reaches the classifier: the same precursors are scored, differently. Whether it helps is
    /// measured on real data (results/entrapment/README.md); the rescorer's own tests show the mechanism.
    /// </summary>
    [Test]
    public void TheClassifierNetworkCanTrainTwice()
    {
        // Decoys elute too, so the classifier has both classes to train on
        var run = SyntheticDiaRun.Build(200, entry => SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(int passes) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(ClassifierNetworkPasses: passes), new CommonParameters(), [], []).Run();

        var once = SearchWith(1);
        var twice = SearchWith(2);

        Assert.That(twice.Matches.Select(m => m.Score), Is.Not.EqualTo(once.Matches.Select(m => m.Score)));
        Assert.That(twice.Matches.Select(m => m.PrecursorIndex), Is.EquivalentTo(once.Matches.Select(m => m.PrecursorIndex)));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DiaLibrarySearchParameters(ClassifierNetworkPasses: 0));
    }

    /// <summary>
    /// The network trains twice by default. With the pooled training sample, two-seed means at a matched entrapment FDP of 1%
    /// gained 2.82% on PXD022589 (HF-X) and 2.00% on PXD005573 1 h (held-out) over one pass: the first network ranks real faint
    /// targets into the pooled sample where the linear model ranked them out.
    /// </summary>
    [Test]
    public void TheClassifierNetworkTrainsTwiceByDefault() =>
        Assert.That(new DiaLibrarySearchParameters().ClassifierNetworkPasses, Is.EqualTo(2));

    /// <summary>
    /// The classifier's ensemble has 24 networks by default, twice DIA-NN's 12. At the 64-32-16 network, two-seed means at a
    /// matched entrapment FDP of 1% gained 0.12% on PXD022589 (HF-X) and 1.60% on PXD005573 1 h (held-out), whose two seeds
    /// gained alike (1.63%, 1.57%). Calibration trains as many, as the A/B ran (DiaIrtSelfCalibrationTests).
    /// </summary>
    [Test]
    public void TheClassifierTrains24NetworksByDefault()
    {
        Assert.That(new DiaLibrarySearchParameters().ClassifierNetworkMembers, Is.EqualTo(24));
    }

    /// <summary>
    /// The tree blend is on by default and never reaches calibration, which only picks anchors. It reorders the final scores:
    /// the same precursors, the same score values, and q-values that still never fall as the score falls.
    /// </summary>
    [Test]
    public void TheTreeBlendIsOnByDefaultAndReordersTheFinalScores()
    {
        Assert.That(new DiaLibrarySearchParameters().ClassifierTreeBlend, Is.True);
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        List<DiaPrecursorMatch> SearchWith(bool blend) => ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(ClassifierNetworkMembers: 3, ClassifierTreeBlend: blend), new CommonParameters(), [], []).Run()).Matches
            .OrderBy(m => m.PrecursorIndex).ToList();

        var plain = SearchWith(false);
        var blended = SearchWith(true);

        Assert.That(plain.Count(m => m.IsDecoy), Is.GreaterThan(150), "enough scored decoys for each fold's trees to learn from");
        Assert.That(blended.Select(m => m.PrecursorIndex), Is.EqualTo(plain.Select(m => m.PrecursorIndex)));
        Assert.That(blended.Select(m => m.Score).Order(), Is.EqualTo(plain.Select(m => m.Score).Order()));
        Assert.That(blended.Select(m => m.Score), Is.Not.EqualTo(plain.Select(m => m.Score)));
        var targets = blended.Where(m => !m.IsDecoy).OrderByDescending(m => m.Score).ToList();
        for (int i = 1; i < targets.Count; i++)
            Assert.That(targets[i].QValue, Is.GreaterThanOrEqualTo(targets[i - 1].QValue), $"at rank {i}");
    }

    /// <summary>
    /// The classifier's network can train on the most confident rows rather than a random sample, as DIA-NN trains after
    /// removing low-confidence identifications. The setting reaches the classifier: with a cap below the rows, the scores change.
    /// </summary>
    [Test]
    public void TheClassifierNetworkCanTrainOnItsMostConfidentRows()
    {
        var run = SyntheticDiaRun.Build(200, entry => SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(StatisticalModels.NetworkTrainingSample sample) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans,
            library, EndpointMap, new DiaLibrarySearchParameters(MaxNetworkTrainingRows: 40, ClassifierNetworkTrainingSample: sample),
            new CommonParameters(), [], []).Run();

        var random = SearchWith(StatisticalModels.NetworkTrainingSample.Random);
        var confident = SearchWith(StatisticalModels.NetworkTrainingSample.Confident);

        Assert.That(confident.Matches.Select(m => m.PrecursorIndex), Is.EquivalentTo(random.Matches.Select(m => m.PrecursorIndex)));
        Assert.That(confident.Matches.Select(m => m.Score), Is.Not.EqualTo(random.Matches.Select(m => m.Score)));
    }

    /// <summary>
    /// By default the network trains on the 120,000 rows per fold the linear model ranks highest, targets and decoys under
    /// one score threshold. Against the half-and-half confident sample at 250,000, two-seed means at a matched entrapment
    /// FDP of 1% gained 3.58% on PXD022589 (HF-X) and 2.23% on PXD005573 1 h (held-out) at 180,000 rows. Under two network
    /// passes 120,000 gained a further 0.95% and 1.06% over 180,000; 90,000 lost 0.44% on HF-X and 250,000 lost on both.
    /// </summary>
    [Test]
    public void TheNetworkTrainsOnItsMostConfidentRowsByDefault()
    {
        var parameters = new DiaLibrarySearchParameters();
        Assert.That(parameters.ClassifierNetworkTrainingSample, Is.EqualTo(StatisticalModels.NetworkTrainingSample.ConfidentPooled));
        Assert.That(parameters.MaxNetworkTrainingRows, Is.EqualTo(120_000));
    }

    /// <summary>
    /// A precursor's other candidate peaks are kept, with their apex and score, so that a miss where another engine chose a
    /// different peak can be told apart: was that peak among our candidates and scored lower, or never a candidate?
    /// </summary>
    [Test]
    public void APrecursorsLosingCandidatePeaksAreKept()
    {
        // Dense noise gives precursors weaker second and third candidate peaks
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, noisePeaksPerScan: 2000);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        // Without the tree blend, which reorders the chosen candidates after the classifier has picked them
        var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(MaxApexCandidates: 3, ClassifierTreeBlend: false), new CommonParameters(), [], []).Run();

        Assert.That(results.LosingCandidates.Values.Sum(v => v.Count), Is.GreaterThan(0));
        foreach (var match in results.Matches)
        {
            if (!results.LosingCandidates.TryGetValue(match.PrecursorIndex, out var others))
                continue;
            Assert.That(others, Has.Count.LessThanOrEqualTo(2));
            Assert.That(others.All(o => o.Score <= match.Score && o.ApexRt.Value != match.ApexRt.Value), $"precursor {match.PrecursorIndex}");
        }
    }

    /// <summary>
    /// Sibling support can come from each other charge state's top candidate peak only. Taken from any of a sibling's
    /// candidates, a decoy gets several chances for a noise peak to sit near its apex; at equal score, decoys on HF-X had
    /// more sibling support than the targets we missed (results/entrapment/README.md, I14).
    /// </summary>
    [Test]
    public void SiblingSupportCanComeFromTheSiblingsTopCandidateOnly()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, noisePeaksPerScan: 2000,
            withCharge3Siblings: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        int sibling = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "SiblingCoElution");
        DiaLibrarySearchResults SearchWith(bool topOnly) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(SiblingTopCandidateOnly: topOnly), new CommonParameters(), [], []).Run();

        var any = SearchWith(false).Matches.ToDictionary(m => m.PrecursorIndex, m => m.Features[sibling]);
        var top = SearchWith(true).Matches.ToDictionary(m => m.PrecursorIndex, m => m.Features[sibling]);

        Assert.That(top.Keys.Intersect(any.Keys).Count(k => top[k] != any[k]), Is.GreaterThan(0));
        // The default since 2026-10-02: HF-X +3.5% at a matched entrapment FDP of 1%, PXD005573 1 h unchanged (-0.3%)
        Assert.That(new DiaLibrarySearchParameters().SiblingTopCandidateOnly, Is.True);
    }

    /// <summary>
    /// DIA-NN 1.8 scores we lacked. MinCorr is co-elution of spike-suppressed traces (each scan the minimum of itself and its
    /// neighbours). NFCorr is the unfragmented precursor m/z in MS2. ShadowCorr is the trace one isotope below each fragment:
    /// is "our" fragment another ion's M+1? All three are 0 unless asked for. A planted precursor's broad peak survives spike
    /// suppression, so its MinCorr is high.
    /// </summary>
    [Test]
    public void DiaNnScoresAreFilledOnlyWhenAskedFor()
    {
        Assert.That(DiaLibrarySearchEngine.SpikeSuppressed([0, 0, 9, 0, 0, 1, 2, 3, 2, 1]),
            Is.EqualTo(new double[] { 0, 0, 0, 0, 0, 0, 1, 2, 1, 1 }));
        int[] scores = new[] { "MinCorr", "NFCorr", "ShadowCorr" }.Select(n => Array.IndexOf(DiaPrecursorMatch.FeatureNames, n)).ToArray();
        Assert.That(scores, Has.All.GreaterThanOrEqualTo(0));

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(DiaNnScores: on), new CommonParameters(), [], []).Run();

        Assert.That(SearchWith(false).Matches.SelectMany(m => scores.Select(i => m.Features[i])), Has.All.EqualTo(0));
        var planted = SearchWith(true).Matches.Where(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence)).ToList();
        Assert.That(planted, Is.Not.Empty);
        Assert.That(planted.Average(m => m.Features[scores[0]]), Is.GreaterThan(0.5));
        // The default since 2026-10-04: with DIA-NN peak finding at 6 candidates, +1.4% HF-X, +3.6% PXD005573 (matched paired FDP 1%)
        Assert.That(new DiaLibrarySearchParameters().DiaNnScores, Is.True);
    }

    /// <summary>
    /// DIA-NN's pSig: each of the six most intense library fragments' share of their summed signal across the peak, in library
    /// rank, so the classifier can set the observed pattern against the library's. 0 unless asked for. On, a planted
    /// precursor's shares sum to 1 and follow the library: the most intense fragment takes more than the sixth.
    /// </summary>
    [Test]
    public void FragmentSignalSharesAreFilledOnlyWhenAskedFor()
    {
        int[] shares = Enumerable.Range(1, 6).Select(k => Array.IndexOf(DiaPrecursorMatch.FeatureNames, $"SignalShare{k}")).ToArray();
        Assert.That(shares, Has.All.GreaterThanOrEqualTo(0));

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(DiaNnSignalShare: on), new CommonParameters(), [], []).Run();

        Assert.That(SearchWith(false).Matches.SelectMany(m => shares.Select(i => m.Features[i])), Has.All.EqualTo(0));
        var planted = SearchWith(true).Matches.Where(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence)).ToList();
        Assert.That(planted, Is.Not.Empty);
        Assert.That(planted.Select(m => shares.Sum(i => m.Features[i])), Has.All.EqualTo(1).Within(1e-9));
        Assert.That(planted.Average(m => m.Features[shares[0]]), Is.GreaterThan(planted.Average(m => m.Features[shares[5]])));
        Assert.That(new DiaLibrarySearchParameters().DiaNnSignalShare, Is.True, "on by default since the three-seed A/B");
    }

    /// <summary>
    /// The linear model that picks each precursor's candidate peak for the network is refit a set number of times (DIA-NN about
    /// 8). 3 by default (a best-single-feature ranking, then two refits); fewer than 1 is refused. More refits still find
    /// what was planted.
    /// </summary>
    [Test]
    public void TheLinearPickCanBeRefitMoreTimes()
    {
        Assert.That(new DiaLibrarySearchParameters().ClassifierLinearIterations, Is.EqualTo(3));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DiaLibrarySearchParameters(ClassifierLinearIterations: 0));

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(ClassifierLinearIterations: 8), new CommonParameters(), [], []).Run();

        int found = results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence));
        Assert.That(found, Is.GreaterThanOrEqualTo((int)Math.Ceiling(0.9 * run.PlantedSequences.Count)));
    }

    /// <summary>
    /// Library fragments beyond the scored top N enter as their own features (DIA-NN's remaining-fragment correlations), so
    /// they add evidence without diluting the core scores: their co-elution with the profile, the share seen at the apex, and
    /// the co-elution weighted by library intensity. 0 unless asked for, and fewer than 0 is refused. A planted precursor's
    /// extra fragments co-elute.
    /// </summary>
    [Test]
    public void FragmentsBeyondTheTopNAreTheirOwnFeatures()
    {
        int[] extra = new[] { "ExtraCoElution", "ExtraMatchedFraction", "ExtraWeightedCoElution" }.Select(n => Array.IndexOf(DiaPrecursorMatch.FeatureNames, n)).ToArray();
        Assert.That(extra, Has.All.GreaterThanOrEqualTo(0));
        // The default since 2026-10-04: fragments 13-24. Two-seed means at a matched paired FDP of 1%, against the same build
        // without them: HF-X 53,859 -> 54,338 (+0.9%), PXD005573 1 h 39,009 -> 39,211 (+0.5%)
        Assert.That(new DiaLibrarySearchParameters().ExtraFragmentCount, Is.EqualTo(12));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DiaLibrarySearchParameters(ExtraFragmentCount: -1));

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(int count) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(TopFragmentCount: 4, ExtraFragmentCount: count), new CommonParameters(), [], []).Run();

        Assert.That(SearchWith(0).Matches.SelectMany(m => extra.Select(i => m.Features[i])), Has.All.EqualTo(0));
        var planted = SearchWith(2).Matches.Where(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence)).ToList();
        Assert.That(planted, Is.Not.Empty);
        Assert.That(planted.Average(m => m.Features[extra[0]]), Is.GreaterThan(0.5));
        Assert.That(planted.Average(m => m.Features[extra[1]]), Is.GreaterThan(0.5));
    }

    /// <summary>
    /// DIA-NN's co-elution per fragment is the best of its correlations at the full tolerance, 0.45x and 0.2x, so a noise peak
    /// at the wide tolerance cannot spoil a fragment whose real peak sits close. As a feature it is 0 unless asked for, and
    /// never below the full-tolerance co-elution of the same fragments.
    /// </summary>
    [Test]
    public void CoElutionCanTakeEachFragmentsBestTolerance()
    {
        int best = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "MaxToleranceCoElution");
        int plain = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "CoElution");
        Assert.That(best, Is.GreaterThanOrEqualTo(0));
        Assert.That(new DiaLibrarySearchParameters().MaxToleranceCoElution, Is.True);

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(MaxToleranceCoElution: on), new CommonParameters(), [], []).Run();

        Assert.That(SearchWith(false).Matches.Select(m => m.Features[best]), Has.All.EqualTo(0));
        var on = SearchWith(true).Matches;
        Assert.That(on.Where(m => m.Features[plain] > 0).Select(m => m.Features[best] - m.Features[plain]), Has.All.GreaterThanOrEqualTo(-1e-12));
        Assert.That(on.Where(m => !m.IsDecoy && run.PlantedSequences.Contains(m.FullSequence)).Average(m => m.Features[best]), Is.GreaterThan(0.5));
    }

    /// <summary>
    /// DIA-NN's fragment rule: a fragment is scored only if it spans at least 3 residues and lies within 200-1800 m/z, so
    /// short ions shared by many peptides (y1, y2, b2) and the crowded low-m/z region do not count as evidence. On by default;
    /// a search still finds what was planted (the synthetic library's 2-residue fragments are left out).
    /// </summary>
    [Test]
    public void FragmentsCanBeLimitedAsDiaNnLimitsThem()
    {
        Omics.SpectralMatch.MslSpectralLibrary.MslFragmentIon Ion(double mz, int number) => new() { Mz = (float)mz, FragmentNumber = number, Intensity = 1 };
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(Ion(450, 3)), Is.True);
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(Ion(450, 2)), Is.False);
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(Ion(199.9, 4)), Is.False);
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(Ion(1800.1, 9)), Is.False);
        // The limits can be set; DIA-NN's are the defaults
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(Ion(250, 3), minimumMz: 300, minimumResidues: 3), Is.False);
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(Ion(450, 3), minimumMz: 200, minimumResidues: 4), Is.False);
        Assert.That(new DiaLibrarySearchParameters().FragmentMinimumMz, Is.EqualTo(200));
        Assert.That(new DiaLibrarySearchParameters().FragmentMinimumResidues, Is.EqualTo(3));
        // A fragment charge limit can be set too (0, the default, sets none)
        var doubly = new Omics.SpectralMatch.MslSpectralLibrary.MslFragmentIon { Mz = 450, FragmentNumber = 5, Charge = 2, Intensity = 1 };
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(doubly), Is.True);
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(doubly, maximumCharge: 1), Is.False);
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(doubly, maximumCharge: 2), Is.True);
        Assert.That(new DiaLibrarySearchParameters().FragmentMaximumCharge, Is.EqualTo(0));
        // Or relative to the precursor: a fragment cannot carry a precursor's whole charge unless its complement is neutral
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(doubly, precursorCharge: 2, chargeBelowPrecursor: true), Is.False);
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(doubly, precursorCharge: 3, chargeBelowPrecursor: true), Is.True);
        Assert.That(new DiaLibrarySearchParameters().FragmentChargeBelowPrecursor, Is.False);
        // Only terminal ions: FragmentNumber is an internal ion's start residue, and a diagnostic ion has no residues (MSL, 010)
        var internalIon = new Omics.SpectralMatch.MslSpectralLibrary.MslFragmentIon { Mz = 450, FragmentNumber = 4, SecondaryFragmentNumber = 9,
            ProductType = Omics.Fragmentation.ProductType.b, SecondaryProductType = Omics.Fragmentation.ProductType.y, Charge = 1, Intensity = 1 };
        var diagnostic = new Omics.SpectralMatch.MslSpectralLibrary.MslFragmentIon { Mz = 450, FragmentNumber = 5, ProductType = Omics.Fragmentation.ProductType.D, Charge = 1, Intensity = 1 };
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(internalIon), Is.False);
        Assert.That(DiaLibrarySearchEngine.IsDiaNnScorable(diagnostic), Is.False);
        // The default since 2026-10-05: two-seed means at a matched paired FDP of 1%, HF-X 54,777 -> 55,218 (+0.8%),
        // PXD005573 1 h 39,052 -> 40,558 (+3.9%)
        Assert.That(new DiaLibrarySearchParameters().DiaNnFragmentFilter, Is.True);

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(DiaNnFragmentFilter: true), new CommonParameters(), [], []).Run();
        int found = results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence));
        Assert.That(found, Is.GreaterThanOrEqualTo((int)Math.Ceiling(0.9 * run.PlantedSequences.Count)));
    }

    /// <summary>
    /// The MS1 tolerance (precursor traces, isotope envelope and MS1 mass error) can be set; 10 ppm by default. A tighter one
    /// keeps random MS1 peaks out of the traces. A search at 10 ppm still finds what was planted, and a precursor's MS1 mass
    /// error never exceeds the tolerance it was read with.
    /// </summary>
    [Test]
    public void TheMs1ToleranceCanBeSet()
    {
        // 5 ppm since 2026-10-06, with the calibrated MS1 offset applied: two-seed means at a matched paired FDP of 1%,
        // +1.6% HF-X, +1.3% PXD005573 against 10 ppm without it (10 since 2026-10-05, 20 before)
        Assert.That(new DiaLibrarySearchParameters().Ms1TolerancePpm, Is.EqualTo(5));
        int ms1Error = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1MassErrorPpm");

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(Ms1TolerancePpm: 10), new CommonParameters(), [], []).Run();

        Assert.That(results.Matches.Select(m => m.Features[ms1Error]), Has.All.LessThanOrEqualTo(10));
        int found = results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence));
        Assert.That(found, Is.GreaterThanOrEqualTo((int)Math.Ceiling(0.9 * run.PlantedSequences.Count)));
    }

    /// <summary>
    /// MS1 co-elution can also be read at half the MS1 tolerance, as DIA-NN scores MS1 at several tolerances: the precursor's
    /// M0 trace keeps only peaks within half the tolerance, correlated with the fragment profile. 0 unless asked for; a
    /// planted precursor's M0 sits at its library m/z, so its tight correlation stays high.
    /// </summary>
    [Test]
    public void Ms1CoElutionCanBeReadAtHalfTheTolerance()
    {
        int tight = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1TightCorrelation");
        Assert.That(tight, Is.GreaterThanOrEqualTo(0));
        Assert.That(new DiaLibrarySearchParameters().Ms1TightCorrelation, Is.True);

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000, withMs1: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(Ms1TightCorrelation: on), new CommonParameters(), [], []).Run();

        Assert.That(SearchWith(false).Matches.Select(m => m.Features[tight]), Has.All.EqualTo(0));
        var planted = SearchWith(true).Matches.Where(m => run.PlantedSequences.Contains(m.FullSequence)).Select(m => m.Features[tight]).Order().ToArray();
        Assert.That(planted, Is.Not.Empty);
        Assert.That(planted[planted.Length / 2], Is.GreaterThan(0.5));
    }

    /// <summary>
    /// Co-elution can also be read on square-root traces, which damp the apex scans so the peak's flanks weigh more: the core
    /// fragments' mean correlation with the square-rooted profile. 0 unless asked for; a planted precursor's fragments share
    /// one elution shape, so it stays high.
    /// </summary>
    [Test]
    public void CoElutionCanBeReadOnSquareRootTraces()
    {
        int sqrt = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "SqrtCoElution");
        Assert.That(sqrt, Is.GreaterThanOrEqualTo(0));
        Assert.That(new DiaLibrarySearchParameters().SqrtCoElution, Is.True);

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(SqrtCoElution: on), new CommonParameters(), [], []).Run();

        Assert.That(SearchWith(false).Matches.Select(m => m.Features[sqrt]), Has.All.EqualTo(0));
        var planted = SearchWith(true).Matches.Where(m => run.PlantedSequences.Contains(m.FullSequence)).Select(m => m.Features[sqrt]).Order().ToArray();
        Assert.That(planted, Is.Not.Empty);
        Assert.That(planted[planted.Length / 2], Is.GreaterThan(0.5));
    }

    /// <summary>
    /// The MS1 envelope read across the peak, not at one scan: M-1 to M3 summed over the co-elution window's MS1 scans,
    /// against the expected pattern with nothing at M-1. A clean envelope scores near 1; a peak one isotope below as large as
    /// M0 (our M0 is then probably another ion's M+1) pulls it well down.
    /// </summary>
    [Test]
    public void PeakEnvelopeCosinePenalisesAPeakOneIsotopeBelow()
    {
        double[] expected = DiaLibrarySearchEngine.ExpectedIsotopes(1500, 4);
        double[] clean = [0, .. expected.Select(e => 1000 * e)];
        double[] shadowed = [1000 * expected[0], .. expected.Select(e => 1000 * e)];

        Assert.That(DiaLibrarySearchEngine.PeakEnvelopeCosine(clean, expected), Is.EqualTo(1).Within(1e-9));
        Assert.That(DiaLibrarySearchEngine.PeakEnvelopeCosine(shadowed, expected), Is.LessThan(0.85));
        Assert.That(DiaLibrarySearchEngine.PeakEnvelopeCosine(new double[5], expected), Is.EqualTo(0));
    }

    /// <summary>
    /// Each match names its apex MS2 scan (dataRepo 003 point 8: a row can point at a real spectrum through a USI): the
    /// scan's retention time is the apex's, and it is an MS2 scan whose isolation window holds the precursor.
    /// </summary>
    [Test]
    public void EachMatchNamesItsApexScan()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, withMs1: true);
        var byNumber = run.Scans.ToDictionary(s => s.OneBasedScanNumber);

        var matches = Search(run).Matches;

        Assert.That(matches, Is.Not.Empty);
        foreach (var match in matches)
        {
            Assert.That(byNumber.TryGetValue(match.ApexScanNumber, out var scan), $"{match.FullSequence}: scan {match.ApexScanNumber}");
            Assert.That(scan!.MsnOrder, Is.EqualTo(2));
            Assert.That(scan.RetentionTime, Is.EqualTo(match.ApexRt.Value));
            Assert.That(scan.IsolationRange.Contains(match.PrecursorMz));
        }
    }

    /// <summary>
    /// The peak envelope as a feature (Ms1PeakEnvelope, on by default): 0 when switched off; on, a planted precursor's M0-M3 follow the
    /// expected pattern on every MS1 scan, so its summed envelope scores near 1.
    /// </summary>
    [Test]
    public void Ms1EnvelopeCanBeReadAcrossThePeak()
    {
        int envelope = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1PeakEnvelopeCosine");
        Assert.That(envelope, Is.GreaterThanOrEqualTo(0));
        Assert.That(new DiaLibrarySearchParameters().Ms1PeakEnvelope, Is.True);

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000, withMs1: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(Ms1PeakEnvelope: on), new CommonParameters(), [], []).Run();

        Assert.That(SearchWith(false).Matches.Select(m => m.Features[envelope]), Has.All.EqualTo(0));
        var planted = SearchWith(true).Matches.Where(m => run.PlantedSequences.Contains(m.FullSequence)).Select(m => m.Features[envelope]).Order().ToArray();
        Assert.That(planted, Is.Not.Empty);
        Assert.That(planted[planted.Length / 2], Is.GreaterThan(0.95));
    }

    /// <summary>
    /// The peak envelope again from tight peaks only (Ms1PeakEnvelopeTight): a second column reading only MS1 peaks within
    /// 0.6x the MS1 tolerance, so the network sees the envelope at two widths while lookups stay at the full one. 0 unless
    /// asked for; on, a planted precursor's isotopes sit at their m/z, so its tight envelope also scores near 1.
    /// </summary>
    [Test]
    public void Ms1EnvelopeCanAlsoBeReadFromTightPeaksOnly()
    {
        int tight = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1PeakEnvelopeTightCosine");
        Assert.That(tight, Is.GreaterThanOrEqualTo(0));
        Assert.That(new DiaLibrarySearchParameters().Ms1PeakEnvelopeTight, Is.True);

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000, withMs1: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on, double ppmOffset = 0) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(Ms1PeakEnvelopeTight: on, Ms1Offset: ppmOffset == 0 ? null : Ms1OffsetModel.Constant(ppmOffset)), new CommonParameters(), [], []).Run();
        double MedianFor(DiaLibrarySearchResults results)
        {
            var values = results.Matches.Where(m => run.PlantedSequences.Contains(m.FullSequence)).Select(m => m.Features[tight]).Order().ToArray();
            Assert.That(values, Is.Not.Empty);
            return values[values.Length / 2];
        }

        Assert.That(SearchWith(false).Matches.Select(m => m.Features[tight]), Has.All.EqualTo(0));
        Assert.That(MedianFor(SearchWith(true)), Is.GreaterThan(0.95));
        Assert.That(MedianFor(SearchWith(true, ppmOffset: 4)), Is.LessThan(0.5), "peaks 4 ppm off are outside 3 ppm, though inside 5");
    }

    /// <summary>
    /// How much precursor there is, not only its shape: log10(1 + M0 summed over the peak's MS1 scans), and that less
    /// log10(1 + the fragments' summed peak area). A real precursor's MS1 scales with its fragments; a match to noise
    /// fragments has little or no MS1 behind it. Off gives zeros; on, planted precursors carry MS1 that unplanted matches lack.
    /// </summary>
    [Test]
    public void Ms1IntensityCanBeAFeature()
    {
        int intensity = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1LogIntensity");
        int ratio = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1ToFragmentLogRatio");
        Assert.That(intensity, Is.GreaterThanOrEqualTo(0));
        Assert.That(ratio, Is.EqualTo(intensity + 1));
        Assert.That(new DiaLibrarySearchParameters().Ms1Intensity, Is.True);

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000, withMs1: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(Ms1Intensity: on), new CommonParameters(), [], []).Run();
        static double Median(IEnumerable<double> values)
        {
            var sorted = values.Order().ToArray();
            Assert.That(sorted, Is.Not.Empty);
            return sorted[sorted.Length / 2];
        }

        var off = SearchWith(false);
        Assert.That(off.Matches.Select(m => m.Features[intensity]), Has.All.EqualTo(0));
        Assert.That(off.Matches.Select(m => m.Features[ratio]), Has.All.EqualTo(0));

        var on = SearchWith(true).Matches;
        var planted = on.Where(m => run.PlantedSequences.Contains(m.FullSequence)).ToList();
        var unplanted = on.Where(m => !run.PlantedSequences.Contains(m.FullSequence)).ToList();
        Assert.That(Median(planted.Select(m => m.Features[intensity])), Is.GreaterThan(Median(unplanted.Select(m => m.Features[intensity])) + 1));
        Assert.That(Median(planted.Select(m => m.Features[ratio])), Is.GreaterThan(Median(unplanted.Select(m => m.Features[ratio]))));
        Assert.That(planted.Select(m => m.Features[intensity] - Math.Log10(1 + m.Quantity) - m.Features[ratio]), Has.All.EqualTo(0).Within(1e-9),
            "the ratio is the MS1 log less the fragments' log area");
    }

    /// <summary>
    /// The MS1 counterpart of PeakSignalFraction: the share of the extraction window's M0 signal inside this peak. A true peak
    /// owns most of its precursor's MS1; an apex on a minor bump does not. Off gives zeros; on, planted precursors score
    /// higher than unplanted matches, and every value is a share.
    /// </summary>
    [Test]
    public void Ms1PeakSignalFractionCanBeAFeature()
    {
        int share = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1PeakSignalFraction");
        Assert.That(share, Is.GreaterThanOrEqualTo(0));
        Assert.That(new DiaLibrarySearchParameters().Ms1PeakSignalFraction, Is.True);

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000, withMs1: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(Ms1PeakSignalFraction: on), new CommonParameters(), [], []).Run();
        static double Median(IEnumerable<double> values)
        {
            var sorted = values.Order().ToArray();
            Assert.That(sorted, Is.Not.Empty);
            return sorted[sorted.Length / 2];
        }

        Assert.That(SearchWith(false).Matches.Select(m => m.Features[share]), Has.All.EqualTo(0));
        var on = SearchWith(true).Matches;
        Assert.That(on.Select(m => m.Features[share]), Has.All.InRange(0, 1));
        double planted = Median(on.Where(m => run.PlantedSequences.Contains(m.FullSequence)).Select(m => m.Features[share]));
        double unplanted = Median(on.Where(m => !run.PlantedSequences.Contains(m.FullSequence)).Select(m => m.Features[share]));
        Assert.That(planted, Is.GreaterThan(0.5));
        Assert.That(planted, Is.GreaterThan(unplanted));
        Assert.That(on.Where(m => !run.PlantedSequences.Contains(m.FullSequence)).Select(m => m.Features[share]), Has.Some.InRange(0.01, 0.9),
            "an unplanted match that sees stray MS1 mostly sees it outside its peak");
    }

    /// <summary>
    /// How many of the peak's MS1 scans show the precursor at all: the share of distinct MS1 scans across the peak with an M0
    /// found. A faint precursor's summed MS1 is small, but it is still there scan after scan; noise is there now and then.
    /// Off gives zeros; on, planted precursors are seen in most MS1 scans of their peak, unplanted matches rarely.
    /// </summary>
    [Test]
    public void Ms1PeakPointsCanBeAFeature()
    {
        int points = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1PeakPointFraction");
        Assert.That(points, Is.GreaterThanOrEqualTo(0));
        Assert.That(new DiaLibrarySearchParameters().Ms1PeakPoints, Is.True);

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000, withMs1: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(Ms1PeakPoints: on), new CommonParameters(), [], []).Run();
        static double Median(IEnumerable<double> values)
        {
            var sorted = values.Order().ToArray();
            Assert.That(sorted, Is.Not.Empty);
            return sorted[sorted.Length / 2];
        }

        Assert.That(SearchWith(false).Matches.Select(m => m.Features[points]), Has.All.EqualTo(0));
        var on = SearchWith(true).Matches;
        Assert.That(on.Select(m => m.Features[points]), Has.All.InRange(0, 1));
        Assert.That(Median(on.Where(m => run.PlantedSequences.Contains(m.FullSequence)).Select(m => m.Features[points])), Is.GreaterThan(0.75), "faint peak edges can miss M0");
        Assert.That(Median(on.Where(m => !run.PlantedSequences.Contains(m.FullSequence)).Select(m => m.Features[points])), Is.LessThan(0.5));
    }

    /// <summary>
    /// A precursor with no fragments beyond the top N has nothing to read as extras; its extra features were 0, the same as
    /// extras that do not co-elute at all, and short peptides DIA-NN finds were scored low for it. With the fallback, such a
    /// precursor's extra features repeat the main fragments' co-elution, matched fraction and weighted co-elution.
    /// </summary>
    [Test]
    public void PrecursorsWithoutExtraFragmentsCanFallBackToTheMainOnes()
    {
        Assert.That(new DiaLibrarySearchParameters().ExtraFragmentFallback, Is.False);
        int Column(string name) => Array.IndexOf(DiaPrecursorMatch.FeatureNames, name);
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        List<DiaPrecursorMatch> SearchWith(bool fallback) => ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(TopFragmentCount: 1000, ExtraFragmentFallback: fallback), new CommonParameters(), [], []).Run()).Matches;

        var off = SearchWith(false);
        Assert.That(off, Is.Not.Empty);
        Assert.That(off.Select(m => m.Features[Column("ExtraCoElution")]), Has.All.EqualTo(0));

        var on = SearchWith(true);
        Assert.That(on.Select(m => m.Features[Column("ExtraCoElution")]), Is.EqualTo(on.Select(m => m.Features[Column("CoElution")])));
        Assert.That(on.Select(m => m.Features[Column("ExtraMatchedFraction")]), Is.EqualTo(on.Select(m => m.Features[Column("MatchedFraction")])));
        Assert.That(on.Select(m => m.Features[Column("ExtraWeightedCoElution")]), Is.EqualTo(on.Select(m => m.Features[Column("WeightedCoElution")])));
        Assert.That(on.Select(m => m.Features[Column("CoElution")]), Has.Some.GreaterThan(0));
    }

    /// <summary>
    /// The MS1 correlations (M0 and M+1 against the fragment profile) can be measured over their own window, as DIA-NN
    /// correlates MS1 over its whole scan window; null keeps the co-elution window. Nothing else changes.
    /// </summary>
    [Test]
    public void TheMs1CorrelationWindowCanBeSet()
    {
        Assert.That(new DiaLibrarySearchParameters().Ms1CorrelationHalfWidthScans, Is.Null);
        int Column(string name) => Array.IndexOf(DiaPrecursorMatch.FeatureNames, name);
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000, withMs1: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        Dictionary<(int, double), double[]> SearchWith(int? window) => ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(Ms1CorrelationHalfWidthScans: window), new CommonParameters(), [], []).Run()).Matches
            .ToDictionary(m => (m.PrecursorIndex, m.ApexRt.Value), m => m.Features);

        var standard = SearchWith(null);
        Assert.That(SearchWith(3), Is.EqualTo(standard));
        // Scores change with the features, so interference removal can keep slightly different matches: compare those in both
        var wide = SearchWith(6);
        var shared = wide.Keys.Where(standard.ContainsKey).ToList();
        Assert.That(shared.Count, Is.GreaterThan(standard.Count * 9 / 10));
        Assert.That(shared.Count(k => wide[k][Column("Ms1Correlation")] != standard[k][Column("Ms1Correlation")]), Is.GreaterThan(0));
        Assert.That(shared.Count(k => wide[k][Column("Ms1IsotopeCorrelation")] != standard[k][Column("Ms1IsotopeCorrelation")]), Is.GreaterThan(0));
        var others = Enumerable.Range(0, DiaPrecursorMatch.FeatureNames.Length)
            .Where(c => c != Column("Ms1Correlation") && c != Column("Ms1IsotopeCorrelation")).ToList();
        var changed = others.Where(c => shared.Any(k => wide[k][c] != standard[k][c])).Select(c => DiaPrecursorMatch.FeatureNames[c]);
        Assert.That(changed, Is.Empty, "only the MS1 correlations change");
    }

    /// <summary>
    /// DIA-NN's pRT: the precursor's library iRT as a feature, so the classifier can weigh the RT error by where on the
    /// gradient it falls. 0 when off.
    /// </summary>
    [Test]
    public void TheLibraryIrtCanBeAFeature()
    {
        Assert.That(new DiaLibrarySearchParameters().LibraryIrtFeature, Is.False);
        int column = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "LibraryIrt");
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        List<DiaPrecursorMatch> SearchWith(bool on) => ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(LibraryIrtFeature: on), new CommonParameters(), [], []).Run()).Matches;

        Assert.That(SearchWith(false).Select(m => m.Features[column]), Has.All.EqualTo(0));
        var on = SearchWith(true);
        Assert.That(on.Select(m => m.Features[column]), Is.EqualTo(on.Select(m => m.LibraryIrt.Value)));
        Assert.That(on.Select(m => m.Features[column]).Distinct().Count(), Is.GreaterThan(1));
    }

    /// <summary>
    /// The co-elution features' window can be set apart from <see cref="DiaLibrarySearchParameters.ApexHalfWidthScans"/>,
    /// which still finds and spaces the candidate apexes (null uses it for both). With the window at 2 and apex finding at 3,
    /// the apexes are the 3/3 search's, and a peak at the same apex in the 2/2 search has the same co-elution.
    /// </summary>
    [Test]
    public void TheCoElutionWindowCanBeSetApartFromApexFinding()
    {
        Assert.That(new DiaLibrarySearchParameters().CoElutionHalfWidthScans, Is.Null);
        int coElution = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "CoElution");
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        Dictionary<(int, double), double> SearchWith(int apex, int? window) => ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(ApexHalfWidthScans: apex, CoElutionHalfWidthScans: window), new CommonParameters(), [], []).Run()).Matches
            .ToDictionary(m => (m.PrecursorIndex, m.ApexRt.Value), m => m.Features[coElution]);

        var wide = SearchWith(3, null);
        var split = SearchWith(3, 2);
        var narrow = SearchWith(2, null);
        Assert.That(split.Keys, Is.EquivalentTo(wide.Keys), "apex finding is unchanged");
        Assert.That(split.Keys.Count(k => split[k] != wide[k]), Is.GreaterThan(0), "the window changes co-elution");
        var shared = split.Keys.Where(narrow.ContainsKey).ToList();
        Assert.That(shared, Is.Not.Empty);
        Assert.That(shared.Select(k => split[k]), Is.EqualTo(shared.Select(k => narrow[k])));
        Assert.That(SearchWith(3, 3), Is.EqualTo(wide));
    }

    /// <summary>
    /// The classifier network's hidden layers can be set; null is 64-32-16, wider than the rescorer's own 25-20-15-10-5
    /// (DIA-NN 2020's). With the pooled sample at 120,000 rows under two passes, two-seed means at a matched entrapment FDP
    /// of 1% gained 1.21% on PXD022589 (HF-X) and 1.09% on PXD005573 1 h (held-out) over 25-20-15-10-5. The default reaches
    /// the classifier, a search with another architecture gives other scores, and a zero-unit layer is refused.
    /// </summary>
    [Test]
    public void TheClassifierNetworkLayersCanBeSet()
    {
        Assert.That(new DiaLibrarySearchParameters().ClassifierNetworkLayers, Is.Null);
        Assert.That(new DiaLibrarySearchParameters().EffectiveClassifierNetworkLayers, Is.EqualTo(new[] { 64, 32, 16 }));
        Assert.That(new DiaLibrarySearchParameters(ClassifierNetworkLayers: [16, 8]).EffectiveClassifierNetworkLayers, Is.EqualTo(new[] { 16, 8 }));

        var run = SyntheticDiaRun.Build(200, entry => SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(int[]? layers) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(ClassifierNetworkLayers: layers, ClassifierNetworkMembers: 3), new CommonParameters(), [], []).Run();

        var standard = SearchWith(null);
        var narrow = SearchWith([16, 8]);

        Assert.That(standard.Matches.Select(m => m.Score), Is.EqualTo(SearchWith([64, 32, 16]).Matches.Select(m => m.Score)));
        Assert.That(narrow.Matches.Select(m => m.Score), Is.Not.EqualTo(standard.Matches.Select(m => m.Score)));
        Assert.Throws<ArgumentOutOfRangeException>(() => SearchWith([16, 0]), "a zero-unit layer is refused, not ignored");
    }

    /// <summary>
    /// The confident network sample can keep only targets that pass a q-value as positives (null, the default, keeps the top
    /// half-cap of targets whatever their q). With a training cap small enough for the confident sample to apply, the scores
    /// change.
    /// </summary>
    [Test]
    public void TheNetworksPositivesCanBeLimitedToPassingTargets()
    {
        Assert.That(new DiaLibrarySearchParameters().ClassifierNetworkPositiveQValue, Is.Null);

        var run = SyntheticDiaRun.Build(200, entry => SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(double? cutoff) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(ClassifierNetworkTrainingSample: StatisticalModels.NetworkTrainingSample.Confident, ClassifierNetworkPositiveQValue: cutoff, MaxNetworkTrainingRows: 100, ClassifierNetworkMembers: 3),
            new CommonParameters(), [], []).Run();

        // Decoys are planted too, so few targets pass a strict cutoff here; a loose one bites
        Assert.That(SearchWith(0.25).Matches.Select(m => m.Score), Is.Not.EqualTo(SearchWith(null).Matches.Select(m => m.Score)));
    }

    /// <summary>
    /// The confident network sample's share of targets can be set (one half by default). With a training cap small enough for
    /// the confident sample to apply, another share gives other scores, and the default share repeats the default search.
    /// </summary>
    [Test]
    public void TheNetworksTargetShareCanBeSet()
    {
        Assert.That(new DiaLibrarySearchParameters().ClassifierNetworkTargetFraction, Is.EqualTo(0.5));

        var run = SyntheticDiaRun.Build(200, entry => SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        List<double> ScoresWith(double? share) => ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            (share is double s ? new DiaLibrarySearchParameters(ClassifierNetworkTrainingSample: StatisticalModels.NetworkTrainingSample.Confident, ClassifierNetworkTargetFraction: s, MaxNetworkTrainingRows: 100, ClassifierNetworkMembers: 3)
                : new DiaLibrarySearchParameters(ClassifierNetworkTrainingSample: StatisticalModels.NetworkTrainingSample.Confident, MaxNetworkTrainingRows: 100, ClassifierNetworkMembers: 3)),
            new CommonParameters(), [], []).Run()).Matches.Select(m => m.Score).ToList();

        var standard = ScoresWith(null);
        Assert.That(ScoresWith(0.5), Is.EqualTo(standard));
        Assert.That(ScoresWith(0.8), Is.Not.EqualTo(standard));
    }

    /// <summary>
    /// MS1 evidence across the peak rather than at the apex: the M0 mass error averaged over the co-elution window, weighted by
    /// intensity, and the M+2 isotope trace's correlation with the fragment profile. Both 0 unless asked for; on, a planted
    /// precursor's M0 sits at its library m/z (small mass error) and its M+2 co-elutes.
    /// </summary>
    [Test]
    public void Ms1EvidenceCanBeReadAcrossThePeak()
    {
        int error = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1PeakMassErrorPpm");
        int m2 = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1Isotope2Correlation");
        Assert.That(new[] { error, m2 }, Has.All.GreaterThanOrEqualTo(0));
        Assert.That(new DiaLibrarySearchParameters().Ms1PeakFeatures, Is.True);

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000, withMs1: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool on) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(Ms1PeakFeatures: on), new CommonParameters(), [], []).Run();

        Assert.That(SearchWith(false).Matches.SelectMany(m => new[] { m.Features[error], m.Features[m2] }), Has.All.EqualTo(0));
        var planted = SearchWith(true).Matches.Where(m => run.PlantedSequences.Contains(m.FullSequence)).ToList();
        double Median(IEnumerable<double> values) { var s = values.Order().ToArray(); return s[s.Length / 2]; }
        Assert.That(planted, Is.Not.Empty);
        Assert.That(Median(planted.Select(m => m.Features[error])), Is.LessThan(1));
        Assert.That(Median(planted.Select(m => m.Features[m2])), Is.GreaterThan(0.5));
    }

    /// <summary>
    /// A run's MS1 offset can be applied: every MS1 lookup moves to the library m/z shifted by the offset at that retention
    /// time, and the MS1 mass error feature is measured from the shifted m/z. Each match also keeps its raw signed apex MS1
    /// error (observed minus library m/z), from which an offset is fitted. On planted precursors, whose M0 sits at the library
    /// m/z, the raw error stays near 0 whatever the offset, and a +3 ppm offset shows up as about 3 ppm of feature error.
    /// </summary>
    [Test]
    public void AnMs1OffsetShiftsTheMs1Lookups()
    {
        Assert.That(new DiaLibrarySearchParameters().Ms1Offset, Is.Null);
        int error = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1MassErrorPpm");
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000, withMs1: true);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        List<DiaPrecursorMatch> PlantedWith(Ms1OffsetModel? offset) => ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library,
            EndpointMap, new DiaLibrarySearchParameters(Ms1Offset: offset), new CommonParameters(), [], []).Run())
            .Matches.Where(m => run.PlantedSequences.Contains(m.FullSequence) && double.IsFinite(m.Ms1ApexErrorPpm)).ToList();
        double Median(IEnumerable<double> values) { var s = values.Order().ToArray(); return s[s.Length / 2]; }

        var none = PlantedWith(null);
        var shifted = PlantedWith(Ms1OffsetModel.Constant(3));

        Assert.That(none, Is.Not.Empty);
        Assert.That(Median(none.Select(m => Math.Abs(m.Ms1ApexErrorPpm))), Is.LessThan(1));
        Assert.That(Median(shifted.Select(m => Math.Abs(m.Ms1ApexErrorPpm))), Is.LessThan(1), "the raw error does not depend on the offset");
        Assert.That(Median(none.Select(m => m.Features[error])), Is.LessThan(1));
        Assert.That(Median(shifted.Select(m => m.Features[error])), Is.EqualTo(3).Within(0.6));
    }

    /// <summary>
    /// The MS1 tolerance can follow the run's own mass accuracy: with Ms1ToleranceSpreadMultiple, it is that multiple of the
    /// calibrated offset's residual spread. Without the option, or without a measured spread, it stays Ms1TolerancePpm.
    /// </summary>
    [Test]
    public void TheMs1ToleranceCanFollowTheRunsMassAccuracy()
    {
        var offset = Ms1OffsetModel.Fit(Enumerable.Range(0, 400).Select(i => (i / 10.0, 2 + (i % 2 == 0 ? 1.0 : -1.0))));
        Assert.That(new DiaLibrarySearchParameters().Ms1ToleranceSpreadMultiple, Is.Null);

        Assert.That(new DiaLibrarySearchParameters(Ms1Offset: offset).EffectiveMs1TolerancePpm, Is.EqualTo(5));
        Assert.That(new DiaLibrarySearchParameters(Ms1Offset: offset, Ms1ToleranceSpreadMultiple: 3).EffectiveMs1TolerancePpm,
            Is.EqualTo(3 * 1.4826).Within(0.1));
        Assert.That(new DiaLibrarySearchParameters(Ms1ToleranceSpreadMultiple: 3).EffectiveMs1TolerancePpm, Is.EqualTo(5), "no offset fitted");
        Assert.That(new DiaLibrarySearchParameters(Ms1Offset: Ms1OffsetModel.Constant(2), Ms1ToleranceSpreadMultiple: 3).EffectiveMs1TolerancePpm,
            Is.EqualTo(5), "no spread measured");
    }

    /// <summary>
    /// A cheap gate before the full features, as DIA-NN's: a scan can be a candidate apex only if at least this many of the
    /// six most intense library fragments are seen there. Noise-only candidates are never scored, so the search scores
    /// fewer of them and still finds what was planted.
    /// </summary>
    [Test]
    public void CandidateApexesCanBeGatedOnTheirTopFragments()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, noisePeaksPerScan: 2000);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(int gate) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(MinimumApexFragments: gate, DiaNnPeakFinding: false), new CommonParameters(), [], []).Run();

        var open = SearchWith(0);
        var gated = SearchWith(4);

        int Candidates(DiaLibrarySearchResults r) => r.Matches.Count + r.LosingCandidates.Values.Sum(v => v.Count);
        Assert.That(Candidates(gated), Is.LessThan(Candidates(open)));
        int found = gated.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence));
        Assert.That(found, Is.GreaterThanOrEqualTo((int)Math.Ceiling(0.9 * run.PlantedSequences.Count)));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DiaLibrarySearchParameters(MinimumApexFragments: 7));
    }

    /// <summary>
    /// One more candidate by a different rule: the apex of the best single fragment's smoothed trace over the whole window,
    /// when no candidate is already there. In about 2,300 of the DIA-NN IDs we miss per file, DIA-NN's peak was never among
    /// our candidates, and more candidates by the same apex score did not help.
    /// </summary>
    [Test]
    public void TheBestFragmentsApexCanBeAnExtraCandidate()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, noisePeaksPerScan: 2000);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(bool extra) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(FragmentApexCandidate: extra), new CommonParameters(), [], []).Run();

        var without = SearchWith(false);
        var with = SearchWith(true);

        int Candidates(DiaLibrarySearchResults r) => r.Matches.Count + r.LosingCandidates.Values.Sum(v => v.Count);
        Assert.That(Candidates(with), Is.GreaterThan(Candidates(without)));
        Assert.That(new DiaLibrarySearchParameters().FragmentApexCandidate, Is.False);
    }

    /// <summary>
    /// The classifier can normalise each fold on a sample of its training groups rather than scoring every training row.
    /// The setting reaches the classifier: the same precursors are scored, and their scores change.
    /// </summary>
    [Test]
    public void TheClassifierCanNormaliseOnASampleOfGroups()
    {
        var run = SyntheticDiaRun.Build(200, entry => SyntheticDiaRun.Bucket(entry, 4) != 0);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults SearchWith(int? groups) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(ClassifierNormalizationGroups: groups), new CommonParameters(), [], []).Run();

        var all = SearchWith(null);
        var sampled = SearchWith(20);

        Assert.That(sampled.Matches.Select(m => m.PrecursorIndex), Is.EquivalentTo(all.Matches.Select(m => m.PrecursorIndex)));
        Assert.That(sampled.Matches.Select(m => m.Score), Is.Not.EqualTo(all.Matches.Select(m => m.Score)));
        // The default since 2026-10-04: each fold's network otherwise scores every row three times over; HF-X 54,017 -> 54,126,
        // PXD005573 1 h 39,282 -> 39,193 at a matched paired FDP of 1%, classifier 3:26 -> 2:51 on HF-X
        Assert.That(new DiaLibrarySearchParameters().ClassifierNormalizationGroups, Is.EqualTo(300_000));
    }

    /// <summary>
    /// A fragment's trace point can be the most intense peak within tolerance, as DIA-NN reads it (diann.cpp 1.8, level()),
    /// rather than the nearest: a faint noise peak closer to the library m/z otherwise stands in for the real fragment.
    /// With none in tolerance there is no peak.
    /// </summary>
    [Test]
    public void AFragmentCanBeReadAsTheMostIntensePeakInTolerance()
    {
        var spectrum = new MassSpectrometry.MzSpectrum([499.990, 500.001, 500.006, 500.03], [5e5, 1e3, 2e5, 9e6], false);
        var tolerance = new MzLibUtil.PpmTolerance(20);

        Assert.That(DiaLibrarySearchEngine.FragmentPeakIndex(spectrum, 500.0, tolerance, mostIntense: false), Is.EqualTo(1));
        Assert.That(DiaLibrarySearchEngine.FragmentPeakIndex(spectrum, 500.0, tolerance, mostIntense: true), Is.EqualTo(0));
        Assert.That(DiaLibrarySearchEngine.FragmentPeakIndex(spectrum, 501.0, tolerance, mostIntense: true), Is.EqualTo(-1));
        Assert.That(DiaLibrarySearchEngine.FragmentPeakIndex(spectrum, 501.0, tolerance, mostIntense: false), Is.EqualTo(-1));
        // The default since 2026-10-04: +0.2% HF-X, +1.4% PXD005573 at a matched paired FDP of 1%
        Assert.That(new DiaLibrarySearchParameters().MostIntenseFragmentPeak, Is.True);
    }

    /// <summary>
    /// DIA-NN 1.8's candidate peaks: a scan is scored by its best fragment's summed correlation to the other top fragments.
    /// It needs at least two fragments, a score of 0.5, and the smoothed best fragment at a local maximum. Every peak within
    /// 2.0 of the best is kept. A five-fragment peak (score near 4) keeps a four-fragment one (near 3) but not a
    /// two-fragment one (near 1). A one-fragment spike is never a candidate, however intense.
    /// </summary>
    [Test]
    public void CandidatePeaksCanBeFoundAsDiaNnFindsThem()
    {
        static double[][] Traces(params (int Apex, int Fragments, double Height)[] peaks)
        {
            var traces = Enumerable.Range(0, 6).Select(_ => new double[40]).ToArray();
            foreach (var (apex, fragments, height) in peaks)
                for (int f = 0; f < fragments; f++)
                    for (int s = Math.Max(0, apex - 5); s <= Math.Min(39, apex + 5); s++)
                        traces[f][s] += height * (f + 1) * Math.Exp(-0.5 * (s - apex) * (s - apex) / 2.25);
            return traces;
        }
        var withFourFragmentPeak = Traces((10, 5, 1e4), (32, 4, 1e4));
        withFourFragmentPeak[0][25] = 1e8;
        var withTwoFragmentPeak = Traces((10, 5, 1e4), (32, 2, 1e4));
        withTwoFragmentPeak[0][25] = 1e8;

        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(withFourFragmentPeak, 3, 10), Is.EqualTo(new[] { 10, 32 }));
        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(withTwoFragmentPeak, 3, 10), Is.EqualTo(new[] { 10 }));
        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(withFourFragmentPeak, 3, 1), Is.EqualTo(new[] { 10 }));
        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(Traces(), 3, 10), Is.Empty);

        // DIA-NN adds the MS1 trace's correlation with the best fragment to the score: two fragments that do not co-elute
        // score 0 and are no candidate, until an MS1 peak shaped like one of them lifts the score to 1
        var gaussian = Enumerable.Range(0, 40).Select(s => Math.Abs(s - 20) <= 5 ? 1e4 * Math.Exp(-0.5 * (s - 20) * (s - 20) / 2.25) : 0).ToArray();
        var inverted = gaussian.Select(v => 2e4 - v).ToArray();
        var unrelated = new[] { gaussian, inverted, new double[40], new double[40], new double[40], new double[40] };
        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(unrelated, 3, 10), Is.Empty);
        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(unrelated, 3, 10, ms1: gaussian), Is.EqualTo(new[] { 20 }));
        Assert.That(new DiaLibrarySearchParameters().DiaNnPeakFindingMs1, Is.False);
        // The default since 2026-10-04, with 6 candidates: neutral at 3, +1.3% HF-X and +2.8% PXD005573 at 6 (matched paired FDP 1%)
        Assert.That(new DiaLibrarySearchParameters().DiaNnPeakFinding, Is.True);
        // 30 since 2026-11-05, with summed-core candidates: two-seed means +0.52% HF-X, +2.43% PXD005573 (matched paired FDP 1%), and targets over the curve gain alike on the
        // species arm (+0.23% / +2.36%, E. coli hits within counting noise). 20 was level with it within noise, and missed the rule on HF-X.
        Assert.That(new DiaLibrarySearchParameters().MaxApexCandidates, Is.EqualTo(30));
    }

    /// <summary>
    /// A real peak the best-fragment rule loses: AADFIDQALAQK 2+ on PXD005573 1 h, at DIA-NN's RT. All six fragments co-elute (score
    /// 4.9 of 5), but the two top scans are close, and each one's best fragment peaks on the other, so neither is a local maximum of its
    /// own best fragment. The sum of the smoothed traces does not depend on which fragment is best, and it peaks at scan 4.
    /// </summary>
    [Test]
    public void APeakWhoseBestFragmentFlipsIsACandidateWhereTheSummedCorePeaks()
    {
        double[][] core =
        [
            [15343, 15918, 35952, 146043, 133197, 41311, 49477, 21132, 0, 0, 0, 0, 0],
            [19104, 36868, 49326, 100855, 92276, 49877, 47894, 40269, 35302, 34846, 14788, 6329, 11153],
            [37112, 44031, 59810, 189770, 190532, 80275, 73164, 66859, 44940, 37505, 23949, 32421, 34708],
            [5068, 3517, 8828, 78810, 66173, 26473, 4181, 0, 0, 0, 0, 0, 0],
            [4437, 6592, 9575, 41817, 41834, 19314, 17656, 11974, 5409, 0, 2700, 0, 0],
            [0, 12479, 12141, 79011, 78581, 30563, 6373, 0, 0, 0, 0, 0, 0],
        ];
        int[] InThePeak(int[] apexes) => apexes.Where(a => a is >= 2 and <= 5).ToArray();

        Assert.That(InThePeak(DiaLibrarySearchEngine.DiaNnCandidateApexes(core, 3, 10)), Is.Empty);
        Assert.That(InThePeak(DiaLibrarySearchEngine.DiaNnCandidateApexes(core, 3, 10, summedCorePeaks: true)), Is.EqualTo(new[] { 4 }));
    }

    /// <summary>
    /// The summed core only adds places a candidate can be: the gates still hold. A one-fragment spike is no candidate though the sum
    /// peaks there, a two-fragment peak stays more than 2.0 below a five-fragment one, and fragments that do not co-elute score below
    /// 0.5 though their sum peaks. Clean peaks are found where they were.
    /// </summary>
    [Test]
    public void TheSummedCoreKeepsTheOtherGatesOfTheCandidateFinder()
    {
        static double[][] Traces(params (int Apex, int Fragments, double Height)[] peaks)
        {
            var traces = Enumerable.Range(0, 6).Select(_ => new double[40]).ToArray();
            foreach (var (apex, fragments, height) in peaks)
                for (int f = 0; f < fragments; f++)
                    for (int s = Math.Max(0, apex - 5); s <= Math.Min(39, apex + 5); s++)
                        traces[f][s] += height * (f + 1) * Math.Exp(-0.5 * (s - apex) * (s - apex) / 2.25);
            return traces;
        }
        var clean = Traces((10, 5, 1e4), (32, 4, 1e4));
        clean[0][25] = 1e8;
        var withTwoFragmentPeak = Traces((10, 5, 1e4), (32, 2, 1e4));
        withTwoFragmentPeak[0][25] = 1e8;
        var gaussian = Enumerable.Range(0, 40).Select(s => Math.Abs(s - 20) <= 5 ? 1e4 * Math.Exp(-0.5 * (s - 20) * (s - 20) / 2.25) : 0).ToArray();
        var opposed = gaussian.Select(v => 2e4 - 0.5 * v).ToArray();
        var sumPeaksButNoCoElution = new[] { gaussian, opposed, new double[40], new double[40], new double[40], new double[40] };

        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(clean, 3, 10, summedCorePeaks: true), Is.EqualTo(new[] { 10, 32 }));
        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(withTwoFragmentPeak, 3, 10, summedCorePeaks: true), Is.EqualTo(new[] { 10 }));
        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(sumPeaksButNoCoElution, 3, 10, summedCorePeaks: true), Is.Empty);

        // A huge interference on one fragment leaves the sum rising with no maximum at the peak. The best-fragment rule still finds it,
        // so the sum adds to that rule and does not replace it.
        var masked = Traces((20, 5, 1e4));
        masked[5] = Enumerable.Range(0, 40).Select(s => 1e7 * s).ToArray();
        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(masked, 3, 10, summedCorePeaks: true), Does.Contain(20));
    }

    /// <summary>
    /// The summed core moves nothing: a sum-peak can only sit more than the half width from every candidate of the best-fragment rule, and
    /// with room for all of them, each of those candidates stays. A first form let a sum-peak with a slightly higher score take over a nearby
    /// apex: on HF-X it moved 3,694 found precursors' apexes by a median 3 scans, and 27% of them were lost, against 2% of those whose apex stayed.
    /// </summary>
    [Test]
    public void TheSummedCoreNeverMovesACandidateOfTheBestFragmentRule()
    {
        var random = new Random(11);
        int added = 0;
        for (int trial = 0; trial < 300; trial++)
        {
            var core = Enumerable.Range(0, 6).Select(_ => new double[60]).ToArray();
            for (int peak = 0; peak < 1 + random.Next(4); peak++)
            {
                int apex = 5 + random.Next(50);
                double height = 1e3 * (1 + random.Next(100)), width = 1 + 2 * random.NextDouble();
                int fragments = 2 + random.Next(5);
                for (int f = 0; f < fragments; f++)
                {
                    double share = 0.2 + random.NextDouble(), jitter = random.NextDouble() - 0.5;
                    for (int s = 0; s < 60; s++)
                        core[f][s] += height * share * Math.Exp(-0.5 * Math.Pow((s - apex - jitter) / width, 2));
                }
            }
            foreach (var trace in core)
                for (int s = 0; s < 60; s++)
                    trace[s] = random.NextDouble() < 0.15 ? 0 : trace[s] * (0.7 + 0.6 * random.NextDouble());

            int[] plain = DiaLibrarySearchEngine.DiaNnCandidateApexes(core, 3, 100);
            int[] summed = DiaLibrarySearchEngine.DiaNnCandidateApexes(core, 3, 100, summedCorePeaks: true);
            int[] capped = DiaLibrarySearchEngine.DiaNnCandidateApexes(core, 3, 2, summedCorePeaks: true);

            Assert.That(summed, Is.SupersetOf(plain), $"trial {trial}");
            foreach (int scan in summed.Concat(capped).Except(plain))
                Assert.That(plain.All(p => Math.Abs(p - scan) > 3), $"trial {trial}: scan {scan} sits on a candidate of the best-fragment rule");
            Assert.That(capped.Length, Is.LessThanOrEqualTo(2), $"trial {trial}");
            added += summed.Length - plain.Length;
        }
        Assert.That(added, Is.GreaterThan(0), "the sum must add candidates for this test to mean anything");
    }

    /// <summary>
    /// A sum-peak competes on score for the capped places. In PXD005573 1 h the peaks this rule rescues mostly belong to precursors whose ten
    /// places the best-fragment rule has already filled, so a sum-peak added after those would never get in, however well it scores. Here
    /// the real flipped peak (score 4.9) beats a weak two-fragment peak (score about 1) for the single place.
    /// </summary>
    [Test]
    public void ASumPeakThatOutscoresACandidateTakesItsPlaceUnderTheCap()
    {
        double[][] real =
        [
            [15343, 15918, 35952, 146043, 133197, 41311, 49477, 21132, 0, 0, 0, 0, 0],
            [19104, 36868, 49326, 100855, 92276, 49877, 47894, 40269, 35302, 34846, 14788, 6329, 11153],
            [37112, 44031, 59810, 189770, 190532, 80275, 73164, 66859, 44940, 37505, 23949, 32421, 34708],
            [5068, 3517, 8828, 78810, 66173, 26473, 4181, 0, 0, 0, 0, 0, 0],
            [4437, 6592, 9575, 41817, 41834, 19314, 17656, 11974, 5409, 0, 2700, 0, 0],
            [0, 12479, 12141, 79011, 78581, 30563, 6373, 0, 0, 0, 0, 0, 0],
        ];
        // Only the peak's own seven scans of the real traces, then a weak two-fragment peak well away from it
        double[][] core = real.Select(trace => trace.Take(8).Concat(new double[32]).ToArray()).ToArray();
        for (int f = 0; f < 2; f++)
            for (int s = 25; s <= 35; s++)
                core[f][s] = 1e4 * (f + 1) * Math.Exp(-0.5 * (s - 30) * (s - 30) / 2.25);

        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(core, 3, 10), Does.Contain(30));
        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(core, 3, 10), Does.Not.Contain(4));
        Assert.That(DiaLibrarySearchEngine.DiaNnCandidateApexes(core, 3, 1, summedCorePeaks: true), Is.EqualTo(new[] { 4 }));
    }

    /// <summary>
    /// The option reaches the finder: on a noisy run it scores more candidate peaks. On by default in the main search since 2026-11-05 (two-seed
    /// means +0.92% HF-X, +1.07% PXD005573 1 h at a matched paired FDP of 1%). Calibration keeps it off, as it was measured.
    /// </summary>
    [Test]
    public void CandidatePeaksFromTheSummedCoreReachTheSearch()
    {
        Assert.That(new DiaLibrarySearchParameters().CandidatePeaksFromSummedCore, Is.True);
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        int Candidates(bool summed)
        {
            var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
                new DiaLibrarySearchParameters(ClassifierNetworkMembers: 3, ClassifierTreeBlend: false, CandidatePeaksFromSummedCore: summed), new CommonParameters(), [], []).Run();
            return results.Matches.Count + results.LosingCandidates.Values.Sum(losers => losers.Count);
        }

        Assert.That(Candidates(true), Is.GreaterThan(Candidates(false)));
    }

    /// <summary>
    /// A fragment's peak can be looked up from a nearby index (the fragment it shadows, one isotope up) instead of by binary
    /// search: the same peak every time, nearest or most intense, from any starting index, including targets off either end.
    /// </summary>
    [Test]
    public void APeakLookedUpFromANearbyIndexIsTheSamePeak()
    {
        var random = new Random(7);
        var tolerance = new MzLibUtil.PpmTolerance(20);
        for (int trial = 0; trial < 200; trial++)
        {
            double[] mz = Enumerable.Range(0, 1 + random.Next(300)).Select(_ => 100 + 1900 * random.NextDouble()).Distinct().OrderBy(x => x).ToArray();
            double[] intensity = mz.Select(_ => (double)random.Next(1, 50)).ToArray();
            var spectrum = new MassSpectrometry.MzSpectrum(mz, intensity, false);
            for (int k = 0; k < 20; k++)
            {
                double target = random.Next(3) == 0 ? mz[random.Next(mz.Length)] * (1 + (random.NextDouble() - 0.5) * 4e-5) : 50 + 2000 * random.NextDouble();
                int hint = random.Next(mz.Length);
                foreach (bool mostIntense in new[] { false, true })
                    Assert.That(DiaLibrarySearchEngine.FragmentPeakIndex(spectrum, target, tolerance, mostIntense, hint),
                        Is.EqualTo(DiaLibrarySearchEngine.FragmentPeakIndex(spectrum, target, tolerance, mostIntense)), $"trial {trial}, target {target}, hint {hint}");
            }
        }
    }

    /// <summary>
    /// A spectrum's m/z bin table: for each 2-Th bin, the first peak at or above the bin's lower edge. It is the start for the
    /// hinted lookup, so a fragment's peak is found without a binary search; any m/z, even off either end, gets a valid start.
    /// </summary>
    [Test]
    public void ASpectrumsBinTableStartsEachBinAtItsFirstPeak()
    {
        var random = new Random(11);
        double[] mz = Enumerable.Range(0, 500).Select(_ => 150 + 1700 * random.NextDouble()).Distinct().OrderBy(x => x).ToArray();
        var spectrum = new MassSpectrometry.MzSpectrum(mz, mz.Select(_ => 1.0).ToArray(), false);
        var bins = DiaLibrarySearchEngine.PeakBins.Of(spectrum);

        for (double edge = 150; edge < 1850; edge += DiaLibrarySearchEngine.PeakBins.Width)
        {
            int start = bins.Hint(edge + 1e-9);
            int first = Array.FindIndex(mz, x => x >= Math.Floor((edge + 1e-9) / DiaLibrarySearchEngine.PeakBins.Width) * DiaLibrarySearchEngine.PeakBins.Width);
            Assert.That(start, Is.EqualTo(first < 0 ? mz.Length - 1 : first), $"bin at {edge}");
        }
        foreach (double target in new[] { 0.0, 100.0, 5000.0 })
            Assert.That(bins.Hint(target), Is.InRange(0, mz.Length - 1));
        Assert.That(DiaLibrarySearchEngine.PeakBins.Of(new MassSpectrometry.MzSpectrum(Array.Empty<double>(), Array.Empty<double>(), false)).Hint(500), Is.EqualTo(0));
    }

    /// <summary>The search reports peptides as well as precursors, one per full sequence, with peptide-level q-values.</summary>
    [Test]
    public void PeptidesAreReportedWithTheirOwnQValues()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);

        var results = Search(run);

        Assert.That(results.Peptides.Select(p => (p.FullSequence, p.IsDecoy)).Distinct().Count(), Is.EqualTo(results.Peptides.Count));
        Assert.That(results.Peptides, Is.EquivalentTo(DiaPeptideFdr.Assign(results.Matches)));
        int found = results.Peptides.Count(p => !p.IsDecoy && p.QValue <= 0.01 && run.PlantedSequences.Contains(p.FullSequence));
        Assert.That(found, Is.GreaterThanOrEqualTo((int)Math.Ceiling(0.9 * run.PlantedSequences.Count)));
        Assert.That(results.ToString(), Does.Contain("Target peptides with q-value <= 0.01:"));
    }

    /// <summary>
    /// Each precursor carries a quantity: its fragments' summed peak areas. It follows abundance, so a precursor planted
    /// at twice the height quantifies at about twice the amount, and one at equal height at about the same.
    /// </summary>
    [Test]
    public void QuantityFollowsAbundance()
    {
        Func<Omics.SpectralMatch.MslSpectralLibrary.MslLibraryEntry, bool> plant = entry => !entry.IsDecoy;
        var baseline = Search(SyntheticDiaRun.Build(100, plant)).Matches.Where(m => !m.IsDecoy && m.QValue <= 0.01)
            .ToDictionary(m => m.PrecursorIndex);
        var doubled = Search(SyntheticDiaRun.Build(100, plant, abundance: e => SyntheticDiaRun.Bucket(e, 2) == 0 ? 2 : 1))
            .Matches.Where(m => !m.IsDecoy && m.QValue <= 0.01).ToList();

        var ratios = doubled.Where(m => baseline.ContainsKey(m.PrecursorIndex))
            .Select(m => (Doubled: SyntheticDiaRun.Bucket(m.FullSequence, 2) == 0, Ratio: m.Quantity / baseline[m.PrecursorIndex].Quantity))
            .ToList();

        Assert.That(baseline.Values.Select(m => m.Quantity), Is.All.GreaterThan(0));
        Assert.That(ratios.Where(r => r.Doubled).Select(r => r.Ratio).ToList(), Has.Count.GreaterThan(20).And.All.InRange(1.8, 2.2));
        Assert.That(ratios.Where(r => !r.Doubled).Select(r => r.Ratio).ToList(), Has.Count.GreaterThan(20).And.All.InRange(0.9, 1.1));
    }

    /// <summary>
    /// Library context as features (DIA-NN feeds m/z, charge and fragment count to its classifier), so the classifier can
    /// learn when to trust a prediction. A decoy shares its target's precursor m/z and charge, so these cannot reveal the label.
    /// </summary>
    [Test]
    public void PrecursorMzChargeAndLibraryFragmentCountAreFeatures()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000,
            withCharge3Siblings: true);
        var entryOf = run.Library.ToDictionary(e => (e.FullSequence, e.ChargeState));

        var matches = Search(run).Matches;

        int mz = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "PrecursorMz");
        int charge = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "PrecursorCharge");
        int count = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "LibraryFragmentCount");
        Assert.That(new[] { mz, charge, count }, Is.All.GreaterThanOrEqualTo(0));
        Assert.That(matches.Select(m => m.Charge).Distinct().Count(), Is.EqualTo(2));
        foreach (var m in matches)
        {
            Assert.That(m.Features[mz], Is.EqualTo(m.PrecursorMz).Within(1e-3));
            Assert.That(m.Features[charge], Is.EqualTo(m.Charge));
            Assert.That(m.Features[count], Is.EqualTo(entryOf[(m.FullSequence, m.Charge)].MatchedFragmentIons.Count));
        }
    }

    /// <summary>
    /// Charge-state siblings (AlphaDIA scores elution groups): a precursor whose other charge state co-elutes at the same
    /// apex has independent evidence. The features use only the sequence and the decoy flag, never which one is a target,
    /// and decoys have sibling pairs too. A precursor with no sibling in the library gets none.
    /// </summary>
    [Test]
    public void ACoElutingChargeStateSiblingIsEvidence()
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000,
            withCharge3Siblings: true);
        var withSibling = run.Library.GroupBy(e => (e.FullSequence, e.IsDecoy)).Where(g => g.Count() > 1).Select(g => g.Key).ToHashSet();

        var matches = Search(run).Matches;

        int sibling = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "SiblingCoElution");
        int delta = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "SiblingApexDeltaMinutes");
        Assert.That(new[] { sibling, delta }, Is.All.GreaterThanOrEqualTo(0));
        double Median(IEnumerable<double> values) { var s = values.Order().ToArray(); return s[s.Length / 2]; }
        var plantedPairs = matches.Where(m => run.PlantedSequences.Contains(m.FullSequence) && withSibling.Contains((m.FullSequence, m.IsDecoy))).ToList();
        var chancePairs = matches.Where(m => !run.PlantedSequences.Contains(m.FullSequence) && withSibling.Contains((m.FullSequence, m.IsDecoy))).ToList();
        var alone = matches.Where(m => !withSibling.Contains((m.FullSequence, m.IsDecoy))).ToList();
        Assert.That(plantedPairs, Has.Count.GreaterThan(30));
        Assert.That(chancePairs, Has.Count.GreaterThan(30));
        Assert.That(chancePairs.Count(m => m.IsDecoy), Is.GreaterThan(10), "decoys have siblings too");
        Assert.That(Median(plantedPairs.Select(m => m.Features[sibling])), Is.GreaterThan(0.8));
        Assert.That(Median(plantedPairs.Select(m => m.Features[delta])), Is.LessThan(0.05));
        Assert.That(Median(chancePairs.Select(m => m.Features[sibling])), Is.LessThan(0.5));
        Assert.That(alone.Select(m => m.Features[sibling]), Is.All.EqualTo(0));
        Assert.That(alone.Select(m => m.Features[delta]), Is.All.EqualTo(1.0), "no sibling counts as a full minute apart");
    }

    /// <summary>
    /// MS1 evidence that does not lean on the fragments. On the whole-proteome library, a fifth of DIA-NN's identifications we
    /// miss show a clean MS1 peak inside DIA-NN's bounds while their fragments are faint, and our only MS1 features were
    /// correlations to the fragment profile. At the apex's MS1 scan: the M0-M3 envelope against the expected isotope pattern,
    /// the M0 mass error, and how much of the window's MS1 trace maximum the apex holds. A decoy has its target's precursor,
    /// so none of these can reveal the label.
    /// </summary>
    [Test]
    public void Ms1EnvelopeMassErrorAndApexShareSeparatePlantedPrecursorsFromChance()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000, withMs1: true);

        var matches = Search(run).Matches;

        int envelope = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1EnvelopeCosine");
        int ppm = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1MassErrorPpm");
        int share = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "Ms1ApexShare");
        Assert.That(new[] { envelope, ppm, share }, Is.All.GreaterThanOrEqualTo(0));
        double Median(IEnumerable<double> values) { var s = values.Order().ToArray(); return s[s.Length / 2]; }
        var planted = matches.Where(m => run.PlantedSequences.Contains(m.FullSequence)).ToList();
        var chance = matches.Where(m => !run.PlantedSequences.Contains(m.FullSequence)).ToList();
        Assert.That(planted, Has.Count.GreaterThan(50));
        Assert.That(chance, Has.Count.GreaterThan(50));
        Assert.That(Median(planted.Select(m => m.Features[envelope])), Is.GreaterThan(0.95));
        Assert.That(Median(chance.Select(m => m.Features[envelope])), Is.LessThan(0.5));
        Assert.That(Median(planted.Select(m => m.Features[ppm])), Is.LessThan(1), "the planted M0 is at the library m/z");
        Assert.That(Median(chance.Select(m => m.Features[ppm])), Is.EqualTo(new DiaLibrarySearchParameters().Ms1TolerancePpm).Within(1e-9), "no M0 found counts as the full tolerance");
        Assert.That(Median(planted.Select(m => m.Features[share])), Is.GreaterThan(0.8));
        Assert.That(matches.Select(m => m.Features[envelope]), Is.All.InRange(0.0, 1.0 + 1e-12));
        Assert.That(matches.Select(m => m.Features[share]), Is.All.InRange(0.0, 1.0));
    }

    /// <summary>
    /// Interference removal drops matches from the reported list, but the results keep them, so a miss can be told
    /// apart from a precursor never scored. Reported and removed together are every precursor a search without removal
    /// reports.
    /// </summary>
    [Test]
    public void MatchesRemovedAsInterferenceAreKeptInTheResults()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000);
        string path = run.WriteLibrary(_directory);
        DiaLibrarySearchResults SearchWith(int explained)
        {
            using var library = MslLibrary.Load(path);
            return (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
                new DiaLibrarySearchParameters(InterferenceExplainedFragments: explained), new CommonParameters(), [], []).Run();
        }

        var with = SearchWith(4);
        var without = SearchWith(0);

        Assert.That(without.RemovedAsInterference, Is.Empty);
        Assert.That(with.Matches.Select(m => m.PrecursorIndex).Intersect(with.RemovedAsInterference.Select(m => m.PrecursorIndex)), Is.Empty);
        Assert.That(with.Matches.Concat(with.RemovedAsInterference).Select(m => m.PrecursorIndex),
            Is.EquivalentTo(without.Matches.Select(m => m.PrecursorIndex)));
    }

    /// <summary>
    /// The main search removes a match once a better co-eluting match explains 3 of its top fragments: against 4, two-seed
    /// means gained 0.75% on HF-X and 1.21% on the held-out run (6 lost 7.4% / 7.7%).
    /// </summary>
    [Test]
    public void InterferenceRemovalNeedsThreeExplainedFragmentsByDefault()
    {
        Assert.That(new DiaLibrarySearchParameters().InterferenceExplainedFragments, Is.EqualTo(3));
    }

    /// <summary>
    /// As DIA-NN, only a better match at the same precursor m/z (or with it as its +1 isotope) can explain another: with the
    /// 3-fragment rule, two-seed means gained 0.61% on HF-X and 0.85% on the held-out run.
    /// </summary>
    [Test]
    public void InterferenceRemovalPairsTheSameMzByDefault()
    {
        Assert.That(new DiaLibrarySearchParameters().InterferenceSameMzOnly, Is.True);
    }

    /// <summary>Dropping a winner whose own ions do not co-elute is off by default, so a search with no floor is the plain rule.</summary>
    [Test]
    public void InterferenceOwnIonFloorIsOffByDefault()
    {
        Assert.That(new DiaLibrarySearchParameters().InterferenceOwnIonFloor, Is.EqualTo(0));
    }

    /// <summary>Anchors saying the run sits 8 iRT above the library: the error scored is measured from there, and the match keeps its library value for the second pass and the anchors.</summary>
    [Test]
    public void TheCompositionModelMovesTheErrorScored()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 300);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var eight = DiaIrtCompositionModel.Fit(run.Library.Where(e => !e.IsDecoy).Select(e => new DiaIrtCompositionModel.Anchor(e.BaseSequence, 8)).ToList());
        int error = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "AbsoluteDeltaIrt");
        List<DiaPrecursorMatch> Planted(DiaIrtCompositionModel? model) => ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(ClassifierNetworkMembers: 3, ClassifierTreeBlend: false, NeighbourEvidence: false, IrtComposition: model), new CommonParameters(), [], []).Run()).Matches
            .Where(m => run.PlantedSequences.Contains(m.FullSequence)).ToList();

        var plain = Planted(null).ToDictionary(m => m.FullSequence);
        var moved = Planted(eight).Where(m => plain.TryGetValue(m.FullSequence, out var p) && Math.Abs(p.ApexRt.Value - m.ApexRt.Value) < 1e-9).ToList();

        Assert.That(moved, Has.Count.GreaterThan(30));
        Assert.That(moved.Select(m => m.Features[error] - plain[m.FullSequence].Features[error]).Order().ElementAt(moved.Count / 2), Is.EqualTo(8).Within(1));
        Assert.That(moved.All(m => m.LibraryIrt.Value == plain[m.FullSequence].LibraryIrt.Value), Is.True, "the match keeps its library value");
    }

    /// <summary>A library misplacing every precursor by 30 iRT: read in a narrow window around its value, nothing is there, and around where the fit expects it, the planted targets are.</summary>
    [Test]
    public void TheCompositionModelMovesTheWindowRead()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 300);
        foreach (var entry in run.Library)
            entry.RetentionTime -= 30;
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var thirty = DiaIrtCompositionModel.Fit(run.Library.Where(e => !e.IsDecoy).Select(e => new DiaIrtCompositionModel.Anchor(e.BaseSequence, 30)).ToList());
        DiaLibrarySearchResults Search(DiaIrtCompositionModel? model) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(IrtHalfWindow: 5, ClassifierNetworkMembers: 3, ClassifierTreeBlend: false, NeighbourEvidence: false, IrtComposition: model), new CommonParameters(), [], []).Run();
        int Found(DiaLibrarySearchResults results) => results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.05 && run.PlantedSequences.Contains(m.FullSequence));

        var moved = Search(thirty);

        Assert.That(Found(Search(null)), Is.LessThan(5));
        Assert.That(Found(moved), Is.GreaterThan(30));
        Assert.That(moved.Matches.Max(m => Math.Abs(m.ApexIrt.Value - m.LibraryIrt.Value - 30)), Is.LessThan(5.5));
    }

    /// <summary>A precursor whose library iRT lies past the run can still be expected inside it, so the search must look past the run for it.</summary>
    [Test]
    public void TheCompositionModelReachesPastTheRun()
    {
        var run = SyntheticDiaRun.Build(400, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 300);
        double top = run.Scans.Where(s => s.MsnOrder == 2).Max(s => EndpointMap.ToIrt(new MzLibUtil.RtMinutes(s.RetentionTime)).Value);
        foreach (var entry in run.Library)
            entry.RetentionTime += 18;
        var past = run.Library.Where(e => !e.IsDecoy && run.PlantedSequences.Contains(e.FullSequence) && e.RetentionTime > top + 5 && e.RetentionTime < top + 15).ToList();
        Assert.That(past, Has.Count.GreaterThan(5), "the case needs planted precursors the query would not reach");
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var back = DiaIrtCompositionModel.Fit(run.Library.Where(e => !e.IsDecoy).Select(e => new DiaIrtCompositionModel.Anchor(e.BaseSequence, -18)).ToList());
        var scored = ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(IrtHalfWindow: 5, ClassifierNetworkMembers: 3, ClassifierTreeBlend: false, NeighbourEvidence: false, IrtComposition: back), new CommonParameters(), [], []).Run())
            .Matches.Select(m => m.FullSequence).ToHashSet();

        Assert.That(past.Count(e => scored.Contains(e.FullSequence)), Is.GreaterThan(past.Count / 2));
    }

    /// <summary>A table that disagrees with the run lowers what a planted target's spectrum reads against the library, so the search does read the library through it.</summary>
    [Test]
    public void TheIntensityTableReachesTheScoring()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 300);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var wrong = DiaFragmentIntensityModel.Fit(run.Library.Where(e => !e.IsDecoy).Select(e =>
            (IReadOnlyList<(MslFragmentIon Ion, double Observed)>)e.MatchedFragmentIons.Select(f => (f, f.Intensity * (f.FragmentNumber % 2 == 0 ? 50.0 : 1.0))).ToList()));
        int cosine = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "ApexCosine");
        double Median(DiaFragmentIntensityModel? table)
        {
            var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
                new DiaLibrarySearchParameters(ClassifierNetworkMembers: 3, ClassifierTreeBlend: false, FragmentIntensities: table), new CommonParameters(), [], []).Run();
            double[] read = results.Matches.Where(m => run.PlantedSequences.Contains(m.FullSequence)).Select(m => m.Features[cosine]).Order().ToArray();
            return read[read.Length / 2];
        }

        Assert.That(Median(wrong), Is.LessThan(Median(null) - 0.01));
    }

    /// <summary>A strong library ion that is not in the run takes a place among the fragments read. Marked down by the table, it gives that place up to one that is there.</summary>
    [Test]
    public void TheIntensityTableChoosesThePrecursorsFragments()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 300);
        foreach (var entry in run.Library)
            entry.MatchedFragmentIons.Add(new MslFragmentIon { Mz = 1450.123f, Intensity = 1, ProductType = Omics.Fragmentation.ProductType.y, FragmentNumber = 14, ResiduePosition = 1, Charge = 1 });
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        var markedDown = DiaFragmentIntensityModel.Fit(run.Library.Where(e => !e.IsDecoy).Select(e =>
            (IReadOnlyList<(MslFragmentIon Ion, double Observed)>)e.MatchedFragmentIons.Select(f => (f, f.FragmentNumber == 14 ? 0.0001 : (double)f.Intensity)).ToList()));
        Assert.That(markedDown.Multiplier(run.Library[0].MatchedFragmentIons[^1]), Is.EqualTo(1 / DiaFragmentIntensityModel.MaximumMultiplier).Within(1e-9));
        int matched = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "MatchedFraction");
        double Mean(DiaFragmentIntensityModel? table) => ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(TopFragmentCount: 3, ExtraFragmentCount: 0, ClassifierNetworkMembers: 3, ClassifierTreeBlend: false, DiaNnFragmentFilter: false, FragmentIntensities: table), new CommonParameters(), [], []).Run())
            .Matches.Where(m => run.PlantedSequences.Contains(m.FullSequence)).Average(m => m.Features[matched]);

        Assert.That(Mean(markedDown), Is.GreaterThan(Mean(null) + 0.05));
    }

    /// <summary>The step sits after the blend and before the q-values, on by default and off in calibration. On, a search still reports the planted targets, and says it ran.</summary>
    [Test]
    public void NeighbourEvidenceRunsInTheSearchByDefault()
    {
        Assert.That(new DiaLibrarySearchParameters().NeighbourEvidence, Is.True);

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 300);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults Search(bool on) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(ClassifierNetworkMembers: 3, ClassifierTreeBlend: false, NeighbourEvidence: on), new CommonParameters(), [], []).Run();
        var lines = new List<string>();
        void Listen(object? sender, StringEventArgs e) => lines.Add(e.S);
        MetaMorpheusEngine.OutLabelStatusHandler += Listen;
        DiaLibrarySearchResults off, on;
        try
        {
            off = Search(false);
            Assert.That(lines, Has.None.Contains("neighbour evidence"));
            on = Search(true);
        }
        finally
        {
            MetaMorpheusEngine.OutLabelStatusHandler -= Listen;
        }

        Assert.That(lines, Has.Some.Contains("neighbour evidence"));
        Assert.That(on.Matches.Where(m => !m.IsDecoy && m.QValue <= 0.01).Select(m => m.FullSequence).Intersect(run.PlantedSequences).Count(),
            Is.GreaterThanOrEqualTo(off.Matches.Where(m => !m.IsDecoy && m.QValue <= 0.01).Select(m => m.FullSequence).Intersect(run.PlantedSequences).Count() * 9 / 10));
    }

    /// <summary>
    /// A run where each of 30 planted targets has a rival on the library: its m/z and iRT and five of its six fragments, with a sixth of its own that is not in the run.
    /// The target's own ion is really there, so it is what tells the two apart. Returns the search at each floor, and the sequences of the target winners with no floor.
    /// </summary>
    private (DiaLibrarySearchResults Plain, DiaLibrarySearchResults Dropping, List<string> TargetWinners) SearchWithRivals(double floor, double shareCeiling = 1)
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000);
        var rivalOf = new Dictionary<string, string>();
        foreach (var target in run.Library.Where(e => !e.IsDecoy && run.PlantedSequences.Contains(e.FullSequence)).Take(30).ToList())
        {
            string sequence = "W" + target.BaseSequence;
            var rival = new MslLibraryEntry { FullSequence = sequence, BaseSequence = sequence, PrecursorMz = target.PrecursorMz, ChargeState = target.ChargeState, RetentionTime = target.RetentionTime };
            foreach (var (fragment, i) in target.MatchedFragmentIons.Select((f, i) => (f, i)))
                rival.MatchedFragmentIons.Add(new MslFragmentIon
                {
                    Mz = i < 5 ? fragment.Mz : fragment.Mz + 3, Intensity = fragment.Intensity, ProductType = fragment.ProductType, FragmentNumber = fragment.FragmentNumber,
                    ResiduePosition = fragment.ResiduePosition, Charge = fragment.Charge,
                });
            run.Library.Add(rival);
            rivalOf[target.FullSequence] = sequence;
        }
        string path = run.WriteLibrary(_directory);
        DiaLibrarySearchResults SearchAt(double f)
        {
            using var library = MslLibrary.Load(path);
            return (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
                new DiaLibrarySearchParameters(InterferenceOwnIonFloor: f, InterferenceOwnIonShareCeiling: shareCeiling), new CommonParameters(), [], []).Run();
        }
        var plain = SearchAt(0);
        var removed = plain.RemovedAsInterference.Select(m => m.FullSequence).ToHashSet();
        var winners = plain.Matches.Where(m => rivalOf.TryGetValue(m.FullSequence, out string? r) && removed.Contains(r)).Select(m => m.FullSequence).ToList();
        return (plain, SearchAt(floor), winners);
    }

    /// <summary>
    /// With the floor on, the engine passes the winners it drops to the results with the rest of the removals: reported and removed together are still every precursor.
    /// A co-elution can never reach a floor of 2, so every winner with a rival goes.
    /// </summary>
    [Test]
    public void WinnersBelowTheOwnIonFloorAreKeptInTheResultsAsRemoved()
    {
        var (plain, dropping, winners) = SearchWithRivals(2);

        Assert.That(winners, Is.Not.Empty, "the case needs a target that beat its rival");
        Assert.That(dropping.Matches.Select(m => m.FullSequence), Does.Not.Contain(winners[0]));
        Assert.That(dropping.RemovedAsInterference.Select(m => m.FullSequence), Does.Contain(winners[0]));
        Assert.That(dropping.Matches.Concat(dropping.RemovedAsInterference).Select(m => m.PrecursorIndex),
            Is.EquivalentTo(plain.Matches.Concat(plain.RemovedAsInterference).Select(m => m.PrecursorIndex)));
    }

    /// <summary>
    /// Only a winner with almost nothing of its own is judged. These targets' sixth fragment carries more than a hundredth of their library intensity, so with the ceiling
    /// below that they are spared, even at a floor no co-elution can reach. With no ceiling the same floor drops them (the test above), so the ceiling is what spares them.
    /// </summary>
    [Test]
    public void AWinnerWithMuchOfItsOwnIsSparedByTheShareCeiling()
    {
        var (_, dropping, winners) = SearchWithRivals(2, shareCeiling: 0.01);

        Assert.That(winners, Is.Not.Empty, "the case needs a target that beat its rival");
        Assert.That(dropping.Matches.Select(m => m.FullSequence), Is.SupersetOf(winners));
    }

    /// <summary>The ceiling spares only what is above it. These targets' one own fragment of six is well under nine tenths of their library intensity, so they are still judged, and they go.</summary>
    [Test]
    public void AWinnerUnderTheShareCeilingIsStillJudged()
    {
        var (_, dropping, winners) = SearchWithRivals(2, shareCeiling: 0.9);

        Assert.That(winners, Is.Not.Empty, "the case needs a target that beat its rival");
        Assert.That(dropping.Matches.Select(m => m.FullSequence), Does.Not.Contain(winners[0]));
    }

    /// <summary>The ceiling is off by default, so the floor alone decides.</summary>
    [Test]
    public void TheOwnIonShareCeilingIsOffByDefault()
    {
        Assert.That(new DiaLibrarySearchParameters().InterferenceOwnIonShareCeiling, Is.EqualTo(1));
    }

    /// <summary>On the real scans a winner's own ion is measured, not assumed: these targets' sixth fragments are in the run on their peak, so none of them goes.</summary>
    [Test]
    public void AWinnerWhoseOwnIonIsInTheRunStays()
    {
        var (_, dropping, winners) = SearchWithRivals(0.1);

        Assert.That(winners, Is.Not.Empty, "the case needs a target that beat its rival");
        Assert.That(dropping.Matches.Select(m => m.FullSequence), Is.SupersetOf(winners));
    }

    /// <summary>
    /// How far apart two apexes may be and still count as co-eluting for interference removal, in cycles: by default the
    /// apex half-width plus one. The removal pass reports the tolerance it used, in minutes.
    /// </summary>
    [Test]
    public void InterferenceRtToleranceCanBeSet()
    {
        Assert.That(new DiaLibrarySearchParameters().InterferenceRtToleranceCycles, Is.Null);
        var run = SyntheticDiaRun.Build(100, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0);
        string path = run.WriteLibrary(_directory);
        double ToleranceMinutes(double? cycles)
        {
            string reported = null;
            EventHandler<StringEventArgs> status = (_, e) => { if (e.S.StartsWith("interference removal:", StringComparison.Ordinal)) reported = e.S; };
            MetaMorpheusEngine.OutLabelStatusHandler += status;
            try
            {
                using var library = MslLibrary.Load(path);
                new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
                    new DiaLibrarySearchParameters(InterferenceRtToleranceCycles: cycles), new CommonParameters(), [], []).Run();
            }
            finally
            {
                MetaMorpheusEngine.OutLabelStatusHandler -= status;
            }
            var match = System.Text.RegularExpressions.Regex.Match(reported, @"RT tolerance ([\d.]+) min");
            return double.Parse(match.Groups[1].Value, System.Globalization.CultureInfo.InvariantCulture);
        }

        double oneCycle = ToleranceMinutes(1);
        Assert.That(oneCycle, Is.GreaterThan(0));
        Assert.That(ToleranceMinutes(null), Is.EqualTo((new DiaLibrarySearchParameters().ApexHalfWidthScans + 1) * oneCycle).Within(2e-4));
        Assert.That(ToleranceMinutes(2.5), Is.EqualTo(2.5 * oneCycle).Within(2e-4));
    }

    /// <summary>
    /// The classifier's network can train on a random subsample of each fold (mzLib's maxNetworkTrainingRows): on a
    /// whole-proteome library, training on every row was two thirds of the search. A cap above the row count changes
    /// nothing; a small one reaches the rescorer.
    /// </summary>
    [Test]
    public void TheNetworkTrainingCapReachesTheClassifier()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000);
        string path = run.WriteLibrary(_directory);
        List<DiaPrecursorMatch> SearchWith(int? cap)
        {
            using var library = MslLibrary.Load(path);
            return ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
                new DiaLibrarySearchParameters(MaxNetworkTrainingRows: cap), new CommonParameters(), [], []).Run()).Matches;
        }

        var all = SearchWith(null);

        Assert.That(SearchWith(10_000_000).Select(m => m.Score), Is.EqualTo(all.Select(m => m.Score)));
        Assert.That(SearchWith(50).Select(m => m.Score), Is.Not.EqualTo(all.Select(m => m.Score)));
    }

    /// <summary>
    /// Peptide length is a feature: a chance match is harder to make with more residues (more fragments, more of them at
    /// sparse high m/z), and on the whole-proteome library the decoys scoring like real identifications were shorter
    /// (median 9 residues against 11). A reversed decoy has exactly its target's length, so the feature cannot reveal the
    /// label.
    /// </summary>
    [Test]
    public void PeptideLengthIsAFeatureAndADecoyHasItsTargetsLength()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000);
        var lengthOf = run.Library.ToDictionary(e => e.FullSequence, e => e.BaseSequence.Length);

        var matches = Search(run).Matches;

        int length = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "PeptideLength");
        Assert.That(length, Is.GreaterThanOrEqualTo(0));
        Assert.That(matches, Has.Count.GreaterThan(100));
        Assert.That(matches.All(m => m.Features[length] == lengthOf[m.FullSequence]));
        Assert.That(matches.Where(m => m.IsDecoy).All(m => m.Features[length] == lengthOf[new string(m.FullSequence.Reverse().ToArray())]));
    }

    /// <summary>
    /// Co-elution is measured against the smoothed best of the six most intense library fragments. A planted precursor's
    /// fragments share one elution profile, so its co-elution is near 1. A candidate seeing only noise has none.
    /// </summary>
    [Test]
    public void CoElutionMeasuresAgreementWithTheBestFragment()
    {
        // Dense noise, so candidates that were never planted still match fragments by chance and get scored
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0, noisePeaksPerScan: 3000);

        var matches = Search(run).Matches;

        int coElution = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "CoElution");
        double Median(IEnumerable<double> values) { var s = values.Order().ToArray(); return s[s.Length / 2]; }
        var planted = matches.Where(m => run.PlantedSequences.Contains(m.FullSequence)).Select(m => m.Features[coElution]).ToList();
        var noise = matches.Where(m => !run.PlantedSequences.Contains(m.FullSequence)).Select(m => m.Features[coElution]).ToList();
        Assert.That(planted, Has.Count.GreaterThan(50));
        Assert.That(noise, Has.Count.GreaterThan(50), "the fixture must score chance matches");
        Assert.That(Median(planted), Is.GreaterThan(0.9));
        Assert.That(Median(noise), Is.LessThan(0.5));
        Assert.That(matches.All(m => m.Features[coElution] is >= 0 and <= 1));
    }

    /// <summary>
    /// The search runs candidates in parallel; the result must not depend on the thread count, match for match, in
    /// order, including every feature, score and q-value.
    /// </summary>
    [Test]
    public void ParallelSearchGivesTheSameResultAsOneThread()
    {
        var run = SyntheticDiaRun.Build(200, entry => SyntheticDiaRun.Bucket(entry, 4) != 0, noisePeaksPerScan: 600);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        DiaLibrarySearchResults Run(int threads) => (DiaLibrarySearchResults)new DiaLibrarySearchEngine(run.Scans, library, EndpointMap,
            new DiaLibrarySearchParameters(), new CommonParameters(maxThreadsToUsePerFile: threads), [], []).Run();

        var one = Run(1).Matches;
        var four = Run(4).Matches;

        Assert.That(four, Has.Count.EqualTo(one.Count));
        for (int i = 0; i < one.Count; i++)
        {
            Assert.That(four[i].PrecursorIndex, Is.EqualTo(one[i].PrecursorIndex), $"row {i}");
            Assert.That(four[i].Score, Is.EqualTo(one[i].Score), $"row {i}");
            Assert.That(four[i].QValue, Is.EqualTo(one[i].QValue), $"row {i}");
            Assert.That(four[i].Features, Is.EqualTo(one[i].Features), $"row {i}");
        }
    }

    [Test]
    public void TargetQValuesNeverDecreaseAsScoreFalls()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0);

        var targets = Search(run).Matches.Where(m => !m.IsDecoy).OrderByDescending(m => m.Score).ToList();

        Assert.That(targets, Is.Not.Empty);
        for (int i = 1; i < targets.Count; i++)
            Assert.That(targets[i].QValue, Is.GreaterThanOrEqualTo(targets[i - 1].QValue), $"at rank {i}");
        Assert.That(targets.All(m => m.QValue is >= 0 and <= 1));
    }

    [Test]
    public void TargetAndDecoyCountsAreReported()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);

        var results = Search(run);

        int targets = results.Matches.Count(m => !m.IsDecoy);
        int decoys = results.Matches.Count(m => m.IsDecoy);
        int passing = results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01);
        Assert.That(results.TargetCount, Is.EqualTo(targets));
        Assert.That(results.DecoyCount, Is.EqualTo(decoys));
        Assert.That(results.ToString(), Does.Contain($"Target precursors: {targets}"));
        Assert.That(results.ToString(), Does.Contain($"Decoy precursors: {decoys}"));
        Assert.That(results.ToString(), Does.Contain($"Target precursors with q-value <= 0.01: {passing}"));
    }

    /// <summary>
    /// iRT and run minutes are different quantities. The types keep them apart, with no conversion in either
    /// direction, not even to double, so the only way across is a map. A search told the identity map (iRT = minutes)
    /// looks in the wrong place and finds almost nothing.
    /// </summary>
    [Test]
    public void IrtAndMinutesDoNotMix()
    {
        foreach (var type in new[] { typeof(Irt), typeof(RtMinutes) })
        {
            var conversions = type.GetMethods(BindingFlags.Public | BindingFlags.Static)
                .Where(method => method.Name is "op_Implicit" or "op_Explicit");
            Assert.That(conversions, Is.Empty, $"{type.Name} must not convert implicitly or explicitly");
        }

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        var identity = IrtCalibration.Line((new RtMinutes(0), new Irt(0)), (new RtMinutes(1), new Irt(1)));

        var results = Search(run, identity, out _);

        int found = results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence));
        Assert.That(found, Is.LessThan(run.PlantedSequences.Count / 2));
    }

    /// <summary>The run-to-iRT map is applied to the run. The library is never rewritten in minutes.</summary>
    [Test]
    public void TheLibraryFileIsUnchangedBySearch()
    {
        var run = SyntheticDiaRun.Build(50, entry => !entry.IsDecoy);
        string path = run.WriteLibrary(_directory);
        byte[] before = SHA256.HashData(File.ReadAllBytes(path));

        using (var library = MslLibrary.Load(path))
            new DiaLibrarySearchEngine(run.Scans, library, EndpointMap, new DiaLibrarySearchParameters(), new CommonParameters(), [], []).Run();

        Assert.That(SHA256.HashData(File.ReadAllBytes(path)), Is.EqualTo(before));
    }

    [Test]
    public void ReportedApexIsInIrtAndInMinutes()
    {
        var run = SyntheticDiaRun.Build(100, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);

        var best = Search(run).Matches.Where(m => !m.IsDecoy && run.PlantedSequences.Contains(m.FullSequence))
            .OrderByDescending(m => m.Score).First();

        var entry = run.Library.Single(e => e.FullSequence == best.FullSequence);
        Assert.That(best.LibraryIrt.Value, Is.EqualTo(entry.RetentionTime).Within(1e-3));
        Assert.That(best.ApexRt.Value, Is.EqualTo(SyntheticDiaRun.TrueRtMinutes(entry.RetentionTime)).Within(SyntheticDiaRun.CycleMinutes));
        Assert.That(best.ApexIrt, Is.EqualTo(EndpointMap.ToIrt(best.ApexRt)));
    }

    [Test]
    public void EachMatchCarriesItsLibraryPrecursor()
    {
        var run = SyntheticDiaRun.Build(100, entry => SyntheticDiaRun.Bucket(entry, 4) != 0);

        var matches = Search(run).Matches;

        Assert.That(matches, Is.Not.Empty);
        foreach (var match in matches)
        {
            var entry = run.Library.Single(e => e.FullSequence == match.FullSequence);
            Assert.That(match.Charge, Is.EqualTo(entry.ChargeState), match.FullSequence);
            Assert.That(match.PrecursorMz, Is.EqualTo(entry.PrecursorMz).Within(1e-3), match.FullSequence);
            Assert.That(match.IsDecoy, Is.EqualTo(entry.IsDecoy), match.FullSequence);
            Assert.That(match.LibraryIrt.Value, Is.EqualTo(entry.RetentionTime).Within(1e-3), match.FullSequence);
        }
    }

    /// <summary>
    /// A precursor isolated by two windows is reported once, and every precursor either window scores alone is scored when
    /// both are searched. The second window covers the first shifted by 0.25 m/z and holds the same scans with every
    /// intensity squared. (The skeleton also checked that the better of the two scores is kept; here each search trains
    /// its own classifier, so scores from separate searches are not on one scale.)
    /// </summary>
    [Test]
    public void APrecursorIsolatedByTwoWindowsIsReportedOnce()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        double firstWindowCentre = SyntheticDiaRun.FirstWindowLowMz + SyntheticDiaRun.WindowWidth / 2;
        int scanNumber = run.Scans.Length + 1;
        var shifted = run.Scans.Where(scan => scan.IsolationMz == firstWindowCentre)
            .Select(scan => SyntheticDiaRun.Scan(
                new MzSpectrum(scan.MassSpectrum.XArray.ToArray(), scan.MassSpectrum.YArray.Select(y => y * y).ToArray(), false),
                scanNumber++, scan.RetentionTime, firstWindowCentre + 0.25))
            .ToList();
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var first = Search(library, run.Scans).Select(m => m.PrecursorIndex).ToHashSet();
        var second = Search(library, shifted).Select(m => m.PrecursorIndex).ToHashSet();
        var both = Search(library, run.Scans.Concat(shifted));

        Assert.That(second.Intersect(first), Is.Not.Empty, "the two windows overlap in the precursors they isolate");
        Assert.That(both.Select(m => m.PrecursorIndex), Is.Unique);
        Assert.That(both.Select(m => m.PrecursorIndex), Is.SupersetOf(first.Union(second)));
    }

    /// <summary>
    /// Only MS2 scans with an isolation range are searched. Adding MS2 scans without an isolation range, holding a perfect
    /// match for a decoy that never eluted over seven cycles, changes no match. The same spectra as MS2 scans isolating the
    /// decoy do change them. (The skeleton also added an MS1 scan claiming an isolation window and an empty MS2 scan inside a window. This
    /// engine reads MS1 scans as MS1 evidence and scores chromatograms, where an empty scan is a point of zero signal, so
    /// both are data here.)
    /// </summary>
    [Test]
    public void ScansThatAreNotDiaWindowsChangeNothing()
    {
        var run = SyntheticDiaRun.Build(100, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        var decoy = run.Library.First(entry => entry.IsDecoy);
        var fragments = decoy.MatchedFragmentIons.OrderBy(f => f.Mz).ToList();
        var perfectMatch = new MzSpectrum(fragments.Select(f => (double)f.Mz).ToArray(), fragments.Select(f => (double)f.Intensity).ToArray(), false);
        double decoyRt = SyntheticDiaRun.TrueRtMinutes(decoy.RetentionTime);
        int scanNumber = run.Scans.Length + 1;
        // A whole elution, rising and falling over seven cycles, since a candidate peak needs co-elution over several scans and an apex
        var elution = Enumerable.Range(-3, 7).Select(k => (Rt: decoyRt + k * SyntheticDiaRun.CycleMinutes,
            Spectrum: new MzSpectrum(perfectMatch.XArray.ToArray(), perfectMatch.YArray.Select(y => y * Math.Exp(-k * k / 4.5)).ToArray(), false))).ToList();
        var extra = elution.Select(e => SyntheticDiaRun.Scan(e.Spectrum, scanNumber++, e.Rt, isolationMz: null)).ToList();
        var sameSpectraAsMs2 = elution.Select(e => SyntheticDiaRun.Scan(e.Spectrum, scanNumber++, e.Rt, decoy.PrecursorMz)).ToList();
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        // A match holds its feature vector as an array, which record equality compares by reference, so compare values
        static List<string> Values(List<DiaPrecursorMatch> matches) =>
            matches.Select(m => $"{m.PrecursorIndex} {m.Score:R} {m.QValue:R} {m.ApexRt.Value:R} {string.Join(",", m.Features.Select(f => f.ToString("R")))}").ToList();
        var plain = Values(Search(library, run.Scans));
        var withExtra = Values(Search(library, run.Scans.Concat(extra)));
        var withMs2 = Values(Search(library, run.Scans.Concat(sameSpectraAsMs2)));

        Assert.That(plain, Is.Not.Empty);
        Assert.That(withExtra, Is.EqualTo(plain));
        Assert.That(withMs2, Is.Not.EqualTo(plain), "searched, a scan in this window would change the matches");
    }

    /// <summary>
    /// An entry with no fragments has nothing to score and is skipped. mzLib's validator refuses such an entry, but
    /// <see cref="MslLibrary.Save"/> does not run it, so a library written elsewhere can hold one.
    /// </summary>
    [Test]
    public void AnEntryWithoutFragmentsIsNotScored()
    {
        var run = SyntheticDiaRun.Build(50, entry => !entry.IsDecoy);
        var template = run.Library.First(entry => !entry.IsDecoy);
        var bare = new MslLibraryEntry
        {
            FullSequence = "PEPTIDEPEPTIDEK",
            BaseSequence = "PEPTIDEPEPTIDEK",
            PrecursorMz = template.PrecursorMz,
            ChargeState = template.ChargeState,
            RetentionTime = template.RetentionTime,
            IsDecoy = false,
        };
        string path = Path.Combine(_directory, $"bare_{Guid.NewGuid():N}.msl");
        MslLibrary.Save(path, run.Library.Append(bare).ToList());
        using var library = MslLibrary.Load(path);

        var sequences = Search(library, run.Scans).Select(m => m.FullSequence).ToList();

        Assert.That(sequences, Does.Contain(template.FullSequence), "its neighbour, at the same m/z and iRT, is still scored");
        Assert.That(sequences, Does.Not.Contain(bare.FullSequence));
    }

    /// <summary>
    /// An entry without a library time (NaN, as mzLib will read a missing one: MSL 021) has no window to search, so it is
    /// not scored, rather than scored against every scan with a NaN library iRT. Its neighbour at the same m/z is.
    /// </summary>
    [Test]
    public void AnEntryWithoutALibraryTimeIsNotScored()
    {
        var run = SyntheticDiaRun.Build(50, entry => !entry.IsDecoy);
        var template = run.Library.First(entry => !entry.IsDecoy);
        var untimed = new MslLibraryEntry
        {
            FullSequence = "NANTIME" + template.FullSequence,
            BaseSequence = "NANTIME" + template.BaseSequence,
            PrecursorMz = template.PrecursorMz,
            ChargeState = template.ChargeState,
            RetentionTime = double.NaN,
            IsDecoy = false,
            MatchedFragmentIons = template.MatchedFragmentIons,
        };
        string path = Path.Combine(_directory, $"untimed_{Guid.NewGuid():N}.msl");
        MslLibrary.Save(path, run.Library.Append(untimed).ToList());
        using var library = MslLibrary.Load(path);

        var matches = Search(library, run.Scans);

        Assert.That(matches.Select(m => m.FullSequence), Does.Contain(template.FullSequence), "its neighbour, at the same m/z, is still scored");
        Assert.That(matches.Select(m => m.FullSequence), Does.Not.Contain(untimed.FullSequence));
        Assert.That(matches, Has.All.Matches<DiaPrecursorMatch>(m => double.IsFinite(m.LibraryIrt.Value)));
    }

    /// <summary>Stop, pressed before the search reaches a window, leaves every window unsearched.</summary>
    [Test]
    [NonParallelizable]
    public void AStoppedSearchSearchesNoWindow()
    {
        var run = SyntheticDiaRun.Build(50, entry => !entry.IsDecoy);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        bool stopLoops = GlobalVariables.StopLoops;
        List<DiaPrecursorMatch> matches;
        try
        {
            GlobalVariables.StopLoops = true;
            matches = Search(library, run.Scans);
        }
        finally
        {
            GlobalVariables.StopLoops = stopLoops;
        }

        Assert.That(matches, Is.Empty);
        Assert.That(Search(library, run.Scans), Is.Not.Empty, "the same search, not stopped, finds precursors");
    }

    [TestCase("scans")]
    [TestCase("library")]
    [TestCase("irtMap")]
    [TestCase("parameters")]
    public void AMissingArgumentIsRefusedByName(string argument)
    {
        var run = SyntheticDiaRun.Build(5, _ => false);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var e = Assert.Throws<ArgumentNullException>(() => new DiaLibrarySearchEngine(
            argument == "scans" ? null! : run.Scans,
            argument == "library" ? null! : library,
            argument == "irtMap" ? null! : EndpointMap,
            argument == "parameters" ? null! : new DiaLibrarySearchParameters(),
            new CommonParameters(), [], []));
        Assert.That(e!.ParamName, Is.EqualTo(argument));
    }
}
