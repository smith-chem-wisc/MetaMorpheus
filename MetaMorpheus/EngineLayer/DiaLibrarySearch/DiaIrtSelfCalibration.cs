#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Chromatography.RetentionTimeCalibration;
using MassSpectrometry;
using MzLibUtil;
using Readers.SpectralLibrary;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>A run's retention-time calibration, and the iRT window a search should use with it.</summary>
/// <param name="AnchorCount">Confident first-pass identifications the calibration was fitted on.</param>
/// <param name="Rounds">Calibration rounds run: the first with the provisional line, each later one with the previous fit.</param>
public sealed record DiaIrtCalibration(IrtCalibrationModel Model, double IrtHalfWindow, int AnchorCount, int Rounds = 1)
{
    /// <summary>
    /// The identifications the calibration was fitted on, as (apex RT, library iRT). A fit can look tight and still be wrong
    /// where its anchors are: on PXD022589 a sample's anchors at 13-27 min, where nothing elutes, bent the curve by up to 39 iRT
    /// and cost 20% of the precursors at an ordinary residual SD.
    /// </summary>
    public IReadOnlyList<(RtMinutes ApexRt, Irt LibraryIrt)> Anchors { get; init; } = [];

    /// <summary>
    /// The run's MS1 offset over retention time, fitted from the last round's anchors' raw MS1 errors; null when calibration
    /// ran a single round. A search applies it through <see cref="DiaLibrarySearchParameters.Ms1Offset"/>.
    /// </summary>
    public Ms1OffsetModel? Ms1Offset { get; init; }

    /// <summary>
    /// How many of the last round's anchors had an M0 error that <see cref="Ms1Offset"/> was fitted from; 0 when calibration ran
    /// a single round. Few against <see cref="AnchorCount"/> means the run has little MS1 signal, or an MS1 error past what
    /// calibration could read, so a task can warn that the offset rests on little.
    /// </summary>
    public int Ms1OffsetAnchors { get; init; }

    /// <summary>
    /// What the last round's anchors say of the library's fragment intensities; null when calibration ran a single round. A search applies it through
    /// <see cref="DiaLibrarySearchParameters.FragmentIntensities"/>.
    /// </summary>
    public DiaFragmentIntensityModel? FragmentIntensities { get; init; }

    /// <summary>
    /// What the last round's anchors say of the library's iRT by sequence composition; null when calibration ran a single round. A search applies it through
    /// <see cref="DiaLibrarySearchParameters.IrtComposition"/>.
    /// </summary>
    public DiaIrtCompositionModel? IrtComposition { get; init; }

    /// <summary>
    /// How many searches the first pass ran before it had anchors to fit, counting every sample and window it tried. A run
    /// that calibrates on its first sample takes one; on a whole-proteome library each costs a search of that sample.
    /// </summary>
    public int FirstPassSearches { get; init; }

    /// <summary>The search's parameters with this calibration applied: its iRT window, its MS1 offset, its fragment intensity table and its iRT composition fit.</summary>
    public DiaLibrarySearchParameters ApplyTo(DiaLibrarySearchParameters parameters) =>
        parameters with { IrtHalfWindow = IrtHalfWindow, Ms1Offset = Ms1Offset, FragmentIntensities = FragmentIntensities, IrtComposition = IrtComposition };
}

/// <summary>
/// Calibrates a DIA run onto its library's iRT scale from the run alone, with nothing borrowed from another search.
/// A first pass searches with a provisional line (the library's iRT range laid over the run's MS2 time span) and a wide
/// window. Its confident targets become anchors for mzLib's robust <see cref="IrtCalibration"/>. The fitting lives in
/// mzLib; this class only orchestrates the passes.
/// </summary>
public static class DiaIrtSelfCalibration
{
    /// <summary>The most confident first-pass targets used as anchors; more add time, not accuracy.</summary>
    public const int MaximumAnchors = 2000;

    /// <summary>
    /// Calibration rounds by default. On PXD022589 (HF-X) 3 rounds give 48,978 precursors at 1% in one search, against
    /// 40,335 for 1 round and 48,417 for a full second search (34 min against 11). On PXD005573 1 h they sit within
    /// seed-to-seed noise (33,676 against a 4-seed mean of 34,334).
    /// </summary>
    public const int DefaultRounds = 3;

    /// <summary>The search window's half-width, in residual SDs of the calibration (at least 5 iRT).</summary>
    public const double DefaultWindowSds = 4;

    /// <summary>Anchors must pass this target-decoy q-value.</summary>
    public const double AnchorQValue = 0.01;

    /// <summary>
    /// A first pass with fewer anchors than this samples twice as many precursors and searches again, until it has them
    /// or has searched the whole library. A sparse library (a whole proteome, of which a run holds a few percent) needs this.
    /// </summary>
    public const int DesiredAnchors = 500;

    /// <summary>The first pass's iRT half-window, as a fraction of the library's central iRT range.</summary>
    public const double FirstPassWindowFraction = 0.25;

    /// <summary>
    /// Targets the first pass aims to search. It finds several times <see cref="MaximumAnchors"/> at 1% from these, so a
    /// larger library is sampled down to about this size.
    /// </summary>
    public const int FirstPassTargetCount = 10_000;

    /// <summary>The first pass's <see cref="DiaLibrarySearchParameters.PrecursorSampleStride"/> for a library of this many targets.</summary>
    public static int FirstPassStride(int targetCount, int firstPassTargetCount = FirstPassTargetCount) =>
        Math.Max(1, (int)Math.Ceiling((double)targetCount / firstPassTargetCount));

    /// <summary>
    /// Calibration's interference rule (<see cref="DiaLibrarySearchParameters.InterferenceExplainedFragments"/>): the main
    /// search's 3 was benchmarked with calibration at 4.
    /// </summary>
    public const int CalibrationInterferenceExplainedFragments = 4;

    /// <summary>
    /// Calibration's classifier network (<see cref="DiaLibrarySearchParameters.ClassifierNetworkLayers"/>): the rescorer's
    /// 25-20-15-10-5, DIA-NN 2020's. The main search's wider default was benchmarked with calibration at this one.
    /// </summary>
    public static IReadOnlyList<int> CalibrationNetworkLayers { get; } = [25, 20, 15, 10, 5];

    /// <summary>
    /// The search's parameters as calibration uses them: every fragment scored, whatever the fragment rule, and MS1 read with
    /// no offset, since calibration is what fits the offset, at <see cref="CalibrationMs1TolerancePpm"/> or the search's own
    /// tolerance if that is wider (a user who knows the instrument is off widens calibration too). With the rule in
    /// calibration too, the held-out run (PXD005573 1 h) lost 2.1% (two-seed means at a matched paired entrapment FDP of 1%).
    /// DIA-NN's signal share stays off as well: it was benchmarked in the main search only, and in calibration it moved the
    /// anchors and cost HF-X 1.5% on seed 0. The network is <see cref="CalibrationNetworkLayers"/>. Calibration only picks
    /// anchors, so it has no use for the tree blend or the neighbour evidence, and it stays on the candidate rule and features it was tuned with, without the summed core.
    /// </summary>
    public static DiaLibrarySearchParameters CalibrationParameters(DiaLibrarySearchParameters parameters) =>
        parameters with { DiaNnFragmentFilter = false, Ms1Offset = null, Ms1TolerancePpm = Math.Max(parameters.Ms1TolerancePpm, CalibrationMs1TolerancePpm), DiaNnSignalShare = false,
            InterferenceExplainedFragments = CalibrationInterferenceExplainedFragments, ClassifierNetworkLayers = [.. CalibrationNetworkLayers],
            ClassifierTreeBlend = false, CandidatePeaksFromSummedCore = false, NeighbourEvidence = false, FragmentIntensities = null, IrtComposition = null,
            MaxToleranceCoElution = false, Ms1TightCorrelation = false, SqrtCoElution = false, Ms1PeakEnvelopeTight = false, Ms1PeakSignalFraction = false, Ms1PeakPoints = false, Ms1PeakFeatures = false,
            Ms1InterpolatedTraces = false, Ms1CorrelationSummedCore = false };

    /// <summary>
    /// MS1 tolerance during calibration, before the run's offset is known: wide enough for an error that drifts up to 3 ppm off
    /// zero (PXD022589). The search's tighter default relies on the offset calibration fits.
    /// </summary>
    public const double CalibrationMs1TolerancePpm = 10;

    /// <summary>
    /// Whether a whole-range search's anchors say the provisional line is misplaced: more than half of them elute outside
    /// the window the first pass searched around the line. A whole-range search can find anchors where the window found none
    /// for another reason, since it searches more of a small sample's precursors, enough to reach 1% q; that alone is no
    /// sign the line is wrong.
    /// </summary>
    internal static bool LineIsMisplaced(IReadOnlyList<(RtMinutes ApexRt, Irt LibraryIrt)> anchors, IrtCalibrationModel provisional, double halfWindow) =>
        anchors.Count > 0 && anchors.Count(a => Math.Abs(provisional.ToIrt(a.ApexRt).Value - a.LibraryIrt.Value) > halfWindow) > anchors.Count / 2.0;

    /// <summary>
    /// Calibrates the run onto the library's iRT scale from the run alone: a first pass with a provisional line and a wide
    /// window, then later rounds, each searching with the previous fit. The last round's anchors also fit the run's MS1
    /// offset, its fragment intensity table and its iRT composition fit (<see cref="DiaIrtCalibration.ApplyTo"/>).
    /// </summary>
    /// <param name="rounds">Calibration rounds, the first pass included (<see cref="DefaultRounds"/>).</param>
    /// <param name="sampleOffset">
    /// Which of each pass's samples to draw (<see cref="DiaLibrarySearchParameters.PrecursorSampleOffset"/>, taken modulo the
    /// pass's stride). 0 by default; others repeat calibration on a disjoint sample.
    /// </param>
    /// <exception cref="ArgumentNullException"><paramref name="scans"/>, <paramref name="library"/> or <paramref name="parameters"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="rounds"/> is below 1, or <paramref name="sampleOffset"/> is negative.</exception>
    /// <exception cref="MetaMorpheusException">
    /// The run has no DIA MS2 scans, the library holds no target precursors with a library time, or the first pass finds too
    /// few confident identifications to calibrate on.
    /// </exception>
    public static DiaIrtCalibration Calibrate(MsDataScan[] scans, MslLibrary library, DiaLibrarySearchParameters parameters,
        CommonParameters commonParameters, IrtCalibrationOptions? options = null, int firstPassTargetCount = FirstPassTargetCount,
        int desiredAnchors = DesiredAnchors, int rounds = DefaultRounds, double windowSds = DefaultWindowSds, int sampleOffset = 0)
    {
        ArgumentNullException.ThrowIfNull(scans);
        ArgumentOutOfRangeException.ThrowIfNegative(sampleOffset);
        if (rounds < 1)
            throw new ArgumentOutOfRangeException(nameof(rounds), rounds, "At least one calibration round is needed.");
        ArgumentNullException.ThrowIfNull(library);
        ArgumentNullException.ThrowIfNull(parameters);
        parameters = CalibrationParameters(parameters);
        options ??= new IrtCalibrationOptions();

        double[] ms2Rts = scans.Where(s => s.MsnOrder == 2 && s.IsolationRange is not null).Select(s => s.RetentionTime).ToArray();
        if (ms2Rts.Length == 0)
            throw new MetaMorpheusException("Could not calibrate retention time: the run has no DIA MS2 scans.");

        // Only targets with a library time set the range: a missing one reads as NaN, which a sort puts below every
        // number, so a few percent of them made the first pass's window NaN
        double[] targetIrts = library.QueryMzWindow(0, float.MaxValue).ToArray()
            .Where(e => e.IsDecoy == 0 && float.IsFinite(e.Irt)).Select(e => (double)e.Irt).Order().ToArray();
        if (targetIrts.Length == 0)
            throw new MetaMorpheusException("Could not calibrate retention time: the library holds no target precursors with a library time.");

        double lowIrt = targetIrts[(int)(0.02 * (targetIrts.Length - 1))];
        double highIrt = targetIrts[(int)(0.98 * (targetIrts.Length - 1))];
        var provisional = IrtCalibration.Line((new RtMinutes(ms2Rts.Min()), new Irt(lowIrt)), (new RtMinutes(ms2Rts.Max()), new Irt(highIrt)));

        // Anchors are picked with the linear discriminant: they only need to be confidently right, and the search's network
        // trained three times over cost most of calibration's time. If the whole library still gives too few, the search's
        // own model gets one try, since on weak data the linear model can find none.
        var linear = parameters with { ClassifierModel = StatisticalModels.RescoreModel.LinearDiscriminant };
        int searches = 0;
        List<(RtMinutes, Irt)> Anchors(DiaLibrarySearchParameters model, int sampleStride, double halfWindow)
        {
            searches++;
            return ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(scans, library, provisional,
                    model with { IrtHalfWindow = halfWindow, PrecursorSampleStride = sampleStride, PrecursorSampleOffset = sampleOffset % sampleStride },
                    commonParameters, [], []).Run()).Matches
                .Where(m => !m.IsDecoy && m.QValue <= AnchorQValue)
                .OrderByDescending(m => m.Score)
                .Take(MaximumAnchors)
                .Select(m => (m.ApexRt, m.LibraryIrt))
                .ToList();
        }

        // The first pass, sampling more of the library until it finds enough anchors or has searched all of it. With mayWiden,
        // a sample where the window finds no anchor at all is searched over the library's whole iRT range too. When that
        // finds enough to fit, its anchors decide (LineIsMisplaced): if most elute outside the window, the line is wrong and the
        // pass continues over the whole range, so a run with dead time stops searching samples with a window that cannot work
        // (PXD071205's 5-min runs: about 38 of calibration's 47 min). If not, the empty sample was only too small to reach 1%
        // q; the wide anchors are set aside, no later sample is widened, and calibration goes on exactly as without this.
        double wholeRange = targetIrts[^1] - targetIrts[0];
        (List<(RtMinutes, Irt)> Anchors, int Stride) FirstPass(double halfWindow, bool mayWiden)
        {
            int stride = FirstPassStride(targetIrts.Length, firstPassTargetCount);
            while (true)
            {
                var anchors = Anchors(linear, stride, halfWindow);
                if (mayWiden && anchors.Count == 0)
                {
                    var wide = Anchors(linear, stride, wholeRange);
                    if (wide.Count >= options.MinimumAnchors)
                    {
                        if (LineIsMisplaced(wide, provisional, halfWindow))
                        {
                            anchors = wide;
                            halfWindow = wholeRange;
                        }
                        mayWiden = false;
                    }
                }
                bool enough = anchors.Count >= Math.Min(desiredAnchors, MaximumAnchors);
                if (!enough && stride == 1 && parameters.ClassifierModel != StatisticalModels.RescoreModel.LinearDiscriminant)
                {
                    var searchModel = Anchors(parameters, stride, halfWindow);
                    if (searchModel.Count > anchors.Count)
                        anchors = searchModel;
                }
                if (enough || stride == 1)
                    return (anchors, stride);
                stride = Math.Max(1, stride / 2);
            }
        }

        var (found, foundStride) = FirstPass(FirstPassWindowFraction * (highIrt - lowIrt), mayWiden: true);
        // The provisional line assumes elution spans the run's MS2 time. A run that spends its first minutes before anything
        // elutes, or a short gradient that elutes only part of a whole-proteome library's iRT range, breaks that, and the line
        // then puts eluting precursors outside the window (PXD071205's 5-min runs: 0.6% inside it, 0 anchors). If no sample's
        // whole-range search found enough either, search once more over the whole range from the first sample, which assumes
        // nothing about where elution is. A run that found enough never reaches this, so its calibration is unchanged.
        if (found.Count < options.MinimumAnchors)
            (found, foundStride) = FirstPass(wholeRange, mayWiden: false);
        int firstPassSearches = searches;
        return LaterRounds(Fit(found, options, windowSds), scans, library, linear, commonParameters, options, foundStride, sampleOffset % foundStride, rounds, windowSds)
            with { FirstPassSearches = firstPassSearches };
    }

    /// <summary>
    /// Iterative calibration, as DIA-NN calibrates in rounds: the same sample searched again with the previous round's fit and
    /// window, so anchors come from the whole gradient rather than where the provisional straight line happened to be right.
    /// A round that cannot fit keeps the previous calibration.
    /// </summary>
    private static DiaIrtCalibration LaterRounds(DiaIrtCalibration calibration, MsDataScan[] scans, MslLibrary library,
        DiaLibrarySearchParameters parameters, CommonParameters commonParameters, IrtCalibrationOptions options, int stride, int offset, int rounds, double windowSds)
    {
        var ms1 = scans.Where(s => s.MsnOrder == 1).OrderBy(s => s.RetentionTime).ToArray();
        double[] ms1Rts = ms1.Select(s => s.RetentionTime).ToArray();
        for (int round = 2; round <= rounds; round++)
        {
            var pass = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(scans, library, calibration.Model,
                parameters with { IrtHalfWindow = calibration.IrtHalfWindow, PrecursorSampleStride = stride, PrecursorSampleOffset = offset }, commonParameters, [], []).Run();
            var confident = pass.Matches
                .Where(m => !m.IsDecoy && m.QValue <= AnchorQValue)
                .OrderByDescending(m => m.Score)
                .Take(MaximumAnchors)
                .ToList();
            var anchors = confident.Select(m => (m.ApexRt, m.LibraryIrt)).ToList();
            if (anchors.Count < options.MinimumAnchors)
                break;
            // The same confident targets give the run's MS1 offset (calibration searches without one, so the errors are raw)
            var fitted = Fit(anchors, options, windowSds);
            var (ms1Offset, ms1OffsetAnchors) = FitMs1Offset(confident, ms1, ms1Rts, parameters.Ms1TolerancePpm);
            calibration = fitted with
            {
                Rounds = round,
                IrtComposition = FitIrtComposition(confident, fitted.Model, library),
                Ms1Offset = ms1Offset,
                Ms1OffsetAnchors = ms1OffsetAnchors,
                FragmentIntensities = FitFragmentIntensities(confident, scans, library, parameters),
            };
        }
        return calibration;
    }

    /// <summary>
    /// The run's MS1 offset from the anchors' M0 errors, and how many anchors had one. Calibration's own search reads MS1 at
    /// <paramref name="tolerancePpm"/> with no offset, so it cannot see an error past that, and it cuts off one side of an
    /// error most of the way there: a run 12 ppm off fitted an offset near 0. So when the anchors' M0s, read from their
    /// nearest MS1 scan at the wider <see cref="DiaMassCalibration.SearchPpm"/>, put the run's median error past half the
    /// tolerance, the offset is fitted from that read instead. A run within half, as every benchmark run so far, keeps the
    /// search's own read and its offset is unchanged.
    /// </summary>
    private static (Ms1OffsetModel Offset, int Anchors) FitMs1Offset(IReadOnlyList<DiaPrecursorMatch> confident, MsDataScan[] ms1, double[] ms1Rts, double tolerancePpm)
    {
        var searched = confident.Select(m => (m.ApexRt.Value, m.Ms1ApexErrorPpm)).ToList();
        if (ms1.Length > 0 && tolerancePpm < DiaMassCalibration.SearchPpm)
        {
            var wide = confident.Select(m => (m.ApexRt.Value, DiaMassCalibration.Ms1ErrorPpm(ms1, ms1Rts, m.ApexRt.Value, m.PrecursorMz))).ToList();
            double[] errors = wide.Select(p => p.Item2).Where(double.IsFinite).Order().ToArray();
            if (errors.Length > 0 && Math.Abs(errors[errors.Length / 2]) > tolerancePpm / 2)
                return (Ms1OffsetModel.Fit(wide), errors.Length);
        }
        return (Ms1OffsetModel.Fit(searched), searched.Count(p => double.IsFinite(p.Item2)));
    }

    private static DiaIrtCompositionModel FitIrtComposition(IReadOnlyList<DiaPrecursorMatch> confident, IrtCalibrationModel model, MslLibrary library) =>
        DiaIrtCompositionModel.Fit(confident.Select(m => (Match: m, Entry: library.GetEntry(m.PrecursorIndex, includeFragments: false))).Where(a => a.Entry is not null)
            .Select(a => new DiaIrtCompositionModel.Anchor(a.Entry!.BaseSequence, model.ToIrt(a.Match.ApexRt).Value - a.Match.LibraryIrt.Value)).ToList());

    private static DiaFragmentIntensityModel FitFragmentIntensities(IReadOnlyList<DiaPrecursorMatch> confident, MsDataScan[] scans, MslLibrary library, DiaLibrarySearchParameters parameters)
    {
        var windows = scans.Where(s => s.MsnOrder == 2 && s.IsolationRange is not null).GroupBy(s => (s.IsolationRange.Minimum, s.IsolationRange.Maximum))
            .Select(w => (w.Key.Minimum, w.Key.Maximum, Scans: w.OrderBy(s => s.RetentionTime).ToArray())).ToArray();
        return DiaFragmentIntensityModel.Fit(confident.Select(m => (Match: m, Entry: library.GetEntry(m.PrecursorIndex))).Where(a => a.Entry is not null).Select(a =>
            (IReadOnlyList<(Omics.SpectralMatch.MslSpectralLibrary.MslFragmentIon Ion, double Observed)>)DiaFragmentIntensityModel.ReadAnchor(a.Entry!, a.Match.ApexRt,
                windows.Where(w => a.Match.PrecursorMz >= w.Minimum && a.Match.PrecursorMz < w.Maximum).Select(w => w.Scans).ToList(), parameters.FragmentTolerancePpm, parameters.MostIntenseFragmentPeak)));
    }

    /// <summary>The most confident main-search targets a <see cref="Refine"/> fit uses.</summary>
    public const int MaximumRefinementAnchors = 20_000;

    /// <summary>
    /// The second pass, as DIA-NN refits after its first search: the run's RT->iRT map refitted on the main search's targets at
    /// <see cref="AnchorQValue"/> (the best <see cref="MaximumRefinementAnchors"/> by score), many more than the first pass's.
    /// Decoys and targets above the cutoff never take part. The refit replaces the curve and its window; the rest of
    /// <paramref name="calibration"/> is kept. Its MS1 offset and intensity table do not depend on the curve. Its iRT
    /// composition fit is read against the curve, so it is refitted against the new one; a calibration without one gets none.
    /// </summary>
    /// <param name="calibration">The calibration the main search ran with.</param>
    /// <param name="library">The library searched, for the composition refit's sequences.</param>
    /// <exception cref="ArgumentNullException">An argument other than <paramref name="options"/> is null.</exception>
    /// <exception cref="MetaMorpheusException">The main search has too few confident targets to refit on.</exception>
    public static DiaIrtCalibration Refine(DiaIrtCalibration calibration, IReadOnlyList<DiaPrecursorMatch> matches, MslLibrary library,
        IrtCalibrationOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(calibration);
        ArgumentNullException.ThrowIfNull(matches);
        ArgumentNullException.ThrowIfNull(library);
        options ??= new IrtCalibrationOptions();
        var confident = matches
            .Where(m => !m.IsDecoy && m.QValue <= AnchorQValue)
            .OrderByDescending(m => m.Score).ThenBy(m => m.PrecursorIndex)
            .Take(MaximumRefinementAnchors)
            .ToList();
        if (confident.Count < options.MinimumAnchors)
            throw new MetaMorpheusException($"Could not refit retention time: the main search found only {confident.Count} confident " +
                $"identifications, and at least {options.MinimumAnchors} are needed.");

        var fitted = Fit(confident.Select(m => (m.ApexRt, m.LibraryIrt)).ToList(), options, DefaultWindowSds);
        return fitted with
        {
            Rounds = calibration.Rounds,
            FirstPassSearches = calibration.FirstPassSearches,
            Ms1Offset = calibration.Ms1Offset,
            Ms1OffsetAnchors = calibration.Ms1OffsetAnchors,
            FragmentIntensities = calibration.FragmentIntensities,
            IrtComposition = calibration.IrtComposition is null ? null : FitIrtComposition(confident, fitted.Model, library),
        };
    }

    private static DiaIrtCalibration Fit(System.Collections.Generic.List<(RtMinutes ApexRt, Irt LibraryIrt)> anchors, IrtCalibrationOptions options, double windowSds)
    {
        if (anchors.Count < options.MinimumAnchors)
            throw new MetaMorpheusException($"Could not calibrate retention time: the first pass found only {anchors.Count} confident " +
                $"identifications, and at least {options.MinimumAnchors} are needed.");

        var model = IrtCalibration.Fit(anchors, options);
        return new DiaIrtCalibration(model, Math.Max(5, windowSds * model.ResidualSd), anchors.Count) { Anchors = anchors };
    }
}
