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
    /// What the last round's anchors say of the library's fragment intensities; null when calibration ran a single round. A search applies it through
    /// <see cref="DiaLibrarySearchParameters.FragmentIntensities"/>.
    /// </summary>
    public DiaFragmentIntensityModel? FragmentIntensities { get; init; }

    /// <summary>
    /// What the last round's anchors say of the library's iRT by sequence composition; null when calibration ran a single round. A search applies it through
    /// <see cref="DiaLibrarySearchParameters.IrtComposition"/>.
    /// </summary>
    public DiaIrtCompositionModel? IrtComposition { get; init; }

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

    /// <param name="sampleOffset">
    /// Which of each pass's samples to draw (<see cref="DiaLibrarySearchParameters.PrecursorSampleOffset"/>, taken modulo the
    /// pass's stride). 0 by default; others repeat calibration on a disjoint sample.
    /// </param>
    /// <exception cref="MetaMorpheusException">
    /// The run has no DIA MS2 scans, the library holds no targets, or too few confident first-pass identifications to
    /// calibrate on.
    /// </exception>
    /// <summary>
    /// The search's parameters as calibration uses them: every fragment scored, whatever the fragment rule, and MS1 read at
    /// <see cref="CalibrationMs1TolerancePpm"/> with no offset, since calibration is what fits the offset. With the rule in
    /// calibration too, the held-out run (PXD005573 1 h) lost 2.1% (two-seed means at a matched paired entrapment FDP of 1%).
    /// DIA-NN's signal share stays off as well: it was benchmarked in the main search only, and in calibration it moved the
    /// anchors and cost HF-X 1.5% on seed 0. The network is <see cref="CalibrationNetworkLayers"/>. Calibration only picks
    /// anchors, so it has no use for the tree blend or the neighbour evidence, and it stays on the candidate rule and features it was tuned with, without the summed core.
    /// </summary>
    public static DiaLibrarySearchParameters CalibrationParameters(DiaLibrarySearchParameters parameters) =>
        parameters with { DiaNnFragmentFilter = false, Ms1Offset = null, Ms1TolerancePpm = CalibrationMs1TolerancePpm, DiaNnSignalShare = false,
            InterferenceExplainedFragments = CalibrationInterferenceExplainedFragments, ClassifierNetworkLayers = [.. CalibrationNetworkLayers],
            ClassifierTreeBlend = false, CandidatePeaksFromSummedCore = false, NeighbourEvidence = false, FragmentIntensities = null, IrtComposition = null,
            MaxToleranceCoElution = false, Ms1TightCorrelation = false, SqrtCoElution = false, Ms1PeakEnvelopeTight = false, Ms1PeakSignalFraction = false, Ms1PeakPoints = false, Ms1PeakFeatures = false };

    /// <summary>
    /// MS1 tolerance during calibration, before the run's offset is known: wide enough for an error that drifts up to 3 ppm off
    /// zero (PXD022589). The search's tighter default relies on the offset calibration fits.
    /// </summary>
    public const double CalibrationMs1TolerancePpm = 10;

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

        // Only targets with a library time set the range: a missing one reads as NaN (MSL 021), which a sort puts below every
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
        List<(RtMinutes, Irt)> Anchors(DiaLibrarySearchParameters model, int sampleStride) =>
            ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(scans, library, provisional,
                model with { IrtHalfWindow = FirstPassWindowFraction * (highIrt - lowIrt), PrecursorSampleStride = sampleStride, PrecursorSampleOffset = sampleOffset % sampleStride },
                commonParameters, [], []).Run()).Matches
            .Where(m => !m.IsDecoy && m.QValue <= AnchorQValue)
            .OrderByDescending(m => m.Score)
            .Take(MaximumAnchors)
            .Select(m => (m.ApexRt, m.LibraryIrt))
            .ToList();

        int stride = FirstPassStride(targetIrts.Length, firstPassTargetCount);
        while (true)
        {
            var anchors = Anchors(linear, stride);
            bool enough = anchors.Count >= Math.Min(desiredAnchors, MaximumAnchors);
            if (!enough && stride == 1 && parameters.ClassifierModel != StatisticalModels.RescoreModel.LinearDiscriminant)
            {
                var searchModel = Anchors(parameters, stride);
                if (searchModel.Count > anchors.Count)
                    anchors = searchModel;
            }
            if (enough || stride == 1)
                return LaterRounds(Fit(anchors, options, windowSds), scans, library, linear, commonParameters, options, stride, sampleOffset % stride, rounds, windowSds);
            stride = Math.Max(1, stride / 2);
        }
    }

    /// <summary>
    /// Iterative calibration, as DIA-NN calibrates in rounds: the same sample searched again with the previous round's fit and
    /// window, so anchors come from the whole gradient rather than where the provisional straight line happened to be right.
    /// A round that cannot fit keeps the previous calibration.
    /// </summary>
    private static DiaIrtCalibration LaterRounds(DiaIrtCalibration calibration, MsDataScan[] scans, MslLibrary library,
        DiaLibrarySearchParameters parameters, CommonParameters commonParameters, IrtCalibrationOptions options, int stride, int offset, int rounds, double windowSds)
    {
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
            calibration = fitted with
            {
                Rounds = round,
                IrtComposition = FitIrtComposition(confident, fitted.Model, library),
                Ms1Offset = Ms1OffsetModel.Fit(confident.Select(m => (m.ApexRt.Value, m.Ms1ApexErrorPpm))),
                FragmentIntensities = FitFragmentIntensities(confident, scans, library, parameters),
            };
        }
        return calibration;
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
    /// Decoys and targets above the cutoff never take part.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="matches"/> is null.</exception>
    /// <exception cref="MetaMorpheusException">Too few confident targets to fit.</exception>
    public static DiaIrtCalibration Refine(IReadOnlyList<DiaPrecursorMatch> matches, IrtCalibrationOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(matches);
        var anchors = matches
            .Where(m => !m.IsDecoy && m.QValue <= AnchorQValue)
            .OrderByDescending(m => m.Score).ThenBy(m => m.PrecursorIndex)
            .Take(MaximumRefinementAnchors)
            .Select(m => (m.ApexRt, m.LibraryIrt))
            .ToList();
        return Fit(anchors, options ?? new IrtCalibrationOptions(), DefaultWindowSds);
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
