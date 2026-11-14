#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MassSpectrometry;
using Chromatography.RetentionTimeCalibration;
using MassSpectrometry.MzSpectra;
using MzLibUtil;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using StatisticalModels;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>
/// Library-based DIA search. It scores every library precursor isolated by each DIA window whose library iRT falls
/// within reach of the run, then assigns target-decoy q-values.
/// </summary>
/// <remarks>
/// For each isolation window, the window's retention-time span is mapped into iRT, and the library is queried there
/// with decoys included. For each candidate, its fragments are read from every scan in the window whose iRT is within
/// <see cref="DiaLibrarySearchParameters.IrtHalfWindow"/> of the library iRT. The apex is the scan with the most summed
/// fragment intensity, and the score is the cosine between the apex intensities and the library's. A precursor isolated
/// by several windows keeps its best score. q-values are (D+1)/T among targets, from mzLib's
/// <see cref="DeconvolutionQValueCalculator.AssignQValues"/>.
/// <para>
/// Everything is compared in iRT. The library is never rewritten in minutes.
/// </para>
/// </remarks>
public class DiaLibrarySearchEngine : MetaMorpheusEngine
{
    private readonly MsDataScan[] _scans;
    private readonly MslLibrary _library;
    private readonly IrtCalibrationModel _irtMap;
    private readonly DiaLibrarySearchParameters _parameters;

    /// <summary>The most intense library fragments that define the elution profile.</summary>
    private const int CoreFragmentCount = 6;

    /// <summary>Tight co-elution counts only peaks within this fraction of the fragment tolerance (DIA-NN uses 0.45).</summary>
    private const double TightToleranceFraction = 0.45;

    /// <summary>MS1 tolerance for the precursor's elution-profile features (DIA-NN chose 17–22 ppm on PXD005573).</summary>
    /// <summary>An MS1 m/z shifted by the run's MS1 offset at <paramref name="rtMinutes"/> (<see cref="DiaLibrarySearchParameters.Ms1Offset"/>).</summary>
    private double Ms1Target(double mz, double rtMinutes) =>
        _parameters.Ms1Offset is { } offset ? mz * (1 + offset.OffsetPpm(rtMinutes) * 1e-6) : mz;

    private PpmTolerance? _ms1Tolerance;
    /// <summary>The MS1 tolerance (<see cref="DiaLibrarySearchParameters.EffectiveMs1TolerancePpm"/>) for precursor traces and envelopes.</summary>
    private PpmTolerance Ms1Tolerance => _ms1Tolerance ??= new PpmTolerance(_parameters.EffectiveMs1TolerancePpm);

    // The run's MS1 peaks, indexed once with FlashLFQ's indexer; null when the run has no MS1 scans
    private FlashLFQ.PeakIndexingEngine? _ms1Index;
    private double[] _ms1Rts = [];

    /// <summary>Each MS2 scan's noise level, sorted for once rather than for every candidate peak that peaks in it.</summary>
    private readonly System.Collections.Concurrent.ConcurrentDictionary<MsDataScan, double> _noiseByScan = new(ReferenceEqualityComparer.Instance);

    /// <summary>A scan's noise: its 25th-percentile intensity.</summary>
    private static double NoiseOf(MsDataScan scan) =>
        scan.MassSpectrum.Size > 0 ? scan.MassSpectrum.YArray.Order().ElementAt(scan.MassSpectrum.Size / 4) : 0;

    /// <param name="scans">The run's scans. Only MS2 scans with an isolation range are searched.</param>
    /// <param name="library">A loaded library holding targets and decoys. The engine does not dispose it.</param>
    /// <param name="irtMap">This run's calibration from minutes to library iRT (mzLib <see cref="IrtCalibrationModel"/>).</param>
    public DiaLibrarySearchEngine(MsDataScan[] scans, MslLibrary library, IrtCalibrationModel irtMap, DiaLibrarySearchParameters parameters,
        CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters,
        List<string> nestedIds)
        : base(commonParameters, fileSpecificParameters, nestedIds)
    {
        _scans = scans ?? throw new ArgumentNullException(nameof(scans));
        _library = library ?? throw new ArgumentNullException(nameof(library));
        _irtMap = irtMap ?? throw new ArgumentNullException(nameof(irtMap));
        _parameters = parameters ?? throw new ArgumentNullException(nameof(parameters));
    }

    /// <exception cref="MetaMorpheusException">The library holds no decoys, so no q-value could be estimated.</exception>
    protected override MetaMorpheusEngineResults RunSpecific()
    {
        if (_library.DecoyCount == 0)
            throw new MetaMorpheusException("The spectral library contains no decoys, so a DIA search cannot estimate its FDR. " +
                "Search a library that includes decoy precursors.");

        Status("Running DIA library search...");
        var tolerance = new PpmTolerance(_parameters.FragmentTolerancePpm);
        var ms1Scans = _scans.Where(scan => scan.MsnOrder == 1).OrderBy(scan => scan.RetentionTime).ToArray();
        _ms1Index = ms1Scans.Length > 0 ? FlashLFQ.PeakIndexingEngine.InitializeIndexingEngine(ms1Scans) : null;
        _ms1Rts = ms1Scans.Select(scan => scan.RetentionTime).ToArray();
        var windows = _scans
            .Where(scan => scan.MsnOrder == 2 && scan.IsolationRange is not null)
            .GroupBy(scan => (scan.IsolationRange.Minimum, scan.IsolationRange.Maximum))
            .OrderBy(window => window.Key.Minimum)
            .ToList();

        // One work item per (window, candidate), scored in parallel; results are gathered in work-item order, so the
        // output does not depend on the thread count
        var work = new List<(MsDataScan[] Scans, double[] ScanIrts, MslPrecursorIndexEntry Candidate)>();
        foreach (var window in windows)
        {
            var scans = window.OrderBy(scan => scan.RetentionTime).ToArray();
            double[] scanIrts = scans.Select(scan => _irtMap.ToIrt(new RtMinutes(scan.RetentionTime)).Value).ToArray();
            double margin = _parameters.IrtComposition is null ? 0 : DiaIrtCompositionModel.QueryMargin;
            double irtLow = scanIrts.Min() - _parameters.IrtHalfWindow - margin;
            double irtHigh = scanIrts.Max() + _parameters.IrtHalfWindow + margin;
            using var hits = _library.QueryWindow((float)window.Key.Minimum, (float)window.Key.Maximum, (float)irtLow, (float)irtHigh, includeDecoys: true);
            foreach (var candidate in hits.Entries.ToArray())
                // A precursor without a library time (NaN, as mzLib reads a missing one) has no window to search
                if (float.IsFinite(candidate.Irt)
                    && candidate.PrecursorIdx % _parameters.PrecursorSampleStride == _parameters.PrecursorSampleOffset % _parameters.PrecursorSampleStride)
                    work.Add((scans, scanIrts, candidate));
        }

        var stage = System.Diagnostics.Stopwatch.StartNew();
        Status($"timing: {work.Count} (window, candidate) work items listed");
        var scored = new List<DiaPrecursorMatch>[work.Count];
        int completed = 0;
        System.Threading.Tasks.Parallel.For(0, work.Count,
            new System.Threading.Tasks.ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUsePerFile },
            (i, loop) =>
            {
                if (GlobalVariables.StopLoops)
                {
                    loop.Stop();
                    return;
                }
                var (scans, scanIrts, candidate) = work[i];
                var entry = _library.GetEntry(candidate.PrecursorIdx);
                scored[i] = entry is null || entry.MatchedFragmentIons.Count == 0
                    ? []
                    : ScoreCandidates(candidate, entry, scans, scanIrts, tolerance).ToList();
                int done = System.Threading.Interlocked.Increment(ref completed);
                if (done % 2000 == 0 || done == work.Count)
                    ReportProgress(new ProgressEventArgs((int)(100.0 * done / work.Count), "Searching DIA candidates...", NestedIds));
            });
        var rows = scored.Where(list => list is not null).SelectMany(list => list).ToList();
        // The lambda above captures these, so they would otherwise live as long as the method: rows alone holds the candidates
        Array.Clear(scored);
        work.Clear();
        work.TrimExcess();
        Status($"timing: extraction and features for {rows.Count} candidate peaks  [{stage.Elapsed:mm\\:ss\\.f}]");

        var windowBounds = windows.Select(w => (Low: w.Key.Minimum, High: w.Key.Maximum)).ToArray();
        double cycleMinutes = windows.Select(w => w.Select(s => s.RetentionTime).Order().ToArray())
            .Where(rts => rts.Length > 1).Select(rts => rts[rts.Length / 2] - rts[rts.Length / 2 - 1]).DefaultIfEmpty(0).Average();
        AddSiblingFeatures(rows, (_parameters.ApexHalfWidthScans + 1) * cycleMinutes, _parameters.SiblingTopCandidateOnly);
        MsDataScan[][] windowScans = _parameters.InterferenceOwnIonFloor > 0 || _parameters.NeighbourEvidence ? windows.Select(w => w.OrderBy(s => s.RetentionTime).ToArray()).ToArray() : [];
        var matches = AssignQValues(rows, windowBounds, windowScans, cycleMinutes, out var removed, out var losing);
        var peptides = DiaPeptideFdr.Assign(matches);
        Status($"timing: search done  [{stage.Elapsed:mm\\:ss\\.f}]");
        Status("Done.");
        return new DiaLibrarySearchResults(this, matches, peptides) { RemovedAsInterference = removed, LosingCandidates = losing };
    }

    /// <summary>
    /// Reads the entry's most intense fragments across the scans within reach of its library iRT, then scores the peak
    /// group with mzLib's <see cref="FragmentCoElution"/>: the apex is where the fragments co-elute in library
    /// proportions, and the score is the apex cosine times the co-elution around it. Up to
    /// <see cref="DiaLibrarySearchParameters.MaxApexCandidates"/> candidate apexes are scored, each as its own match; the
    /// classifier's score decides among them (<see cref="AssignQValues"/>).
    /// </summary>
    private IEnumerable<DiaPrecursorMatch> ScoreCandidates(MslPrecursorIndexEntry candidate, MslLibraryEntry entry, MsDataScan[] scans,
        double[] scanIrts, PpmTolerance tolerance)
    {
        double expectedIrt = candidate.Irt + (_parameters.IrtComposition?.Offset(entry.BaseSequence) ?? 0);
        var allIntensities = entry.MatchedFragmentIons.Select(LibraryIntensity).ToArray();
        // DIA-NN's fragment rule (DiaNnFragmentFilter): ineligible fragments are never picked, unless fewer than 3 qualify
        int eligible = allIntensities.Length;
        if (_parameters.DiaNnFragmentFilter)
        {
            var selection = (double[])allIntensities.Clone();
            int count = 0;
            for (int i = 0; i < selection.Length; i++)
                if (IsDiaNnScorable(entry.MatchedFragmentIons[i], _parameters.FragmentMinimumMz, _parameters.FragmentMinimumResidues, _parameters.FragmentMaximumCharge, entry.ChargeState, _parameters.FragmentChargeBelowPrecursor)) count++;
                else selection[i] = double.NegativeInfinity;
            if (count >= 3)
            {
                allIntensities = selection;
                eligible = count;
            }
        }
        var fragments = FragmentCoElution.TopIndices(allIntensities, Math.Min(_parameters.TopFragmentCount, eligible))
            .Select(i => entry.MatchedFragmentIons[i]).ToList();
        double[] libraryIntensities = fragments.Select(LibraryIntensity).ToArray();
        // Fragments after the top N, by library intensity, read only as their own features (ExtraFragmentCount)
        var extraFragments = _parameters.ExtraFragmentCount == 0 ? [] : Enumerable.Range(0, allIntensities.Length)
            .Where(i => !double.IsNegativeInfinity(allIntensities[i]))
            .OrderByDescending(i => allIntensities[i]).ThenBy(i => i).Skip(fragments.Count).Take(_parameters.ExtraFragmentCount)
            .Select(i => entry.MatchedFragmentIons[i]).ToList();

        // Plain loops in this per-candidate section: it runs for millions of (window, candidate) items, and LINQ's
        // per-element delegates and iterators were a large share of extraction's allocations
        var reachableList = new List<int>();
        for (int s = 0; s < scans.Length; s++)
            if (Math.Abs(scanIrts[s] - expectedIrt) <= _parameters.IrtHalfWindow)
                reachableList.Add(s);
        int[] reachable = reachableList.ToArray();
        if (reachable.Length == 0)
            yield break;

        var traces = new double[fragments.Count][];
        var ppm = new double[fragments.Count][];
        for (int f = 0; f < fragments.Count; f++)
        {
            traces[f] = new double[reachable.Length];
            ppm[f] = new double[reachable.Length];
        }
        // Where each fragment was found, so its isotope shadow (DiaNnScores) is looked up from there, not by binary search
        int[][]? peakIndices = _parameters.DiaNnScores ? fragments.Select(_ => new int[reachable.Length]).ToArray() : null;
        for (int k = 0; k < reachable.Length; k++)
        {
            var spectrum = scans[reachable[k]].MassSpectrum;
            var bins = _binsByScan.GetOrAdd(scans[reachable[k]], s => PeakBins.Of(s.MassSpectrum));
            if (spectrum.Size == 0)
                continue;
            for (int f = 0; f < fragments.Count; f++)
            {
                int i = FragmentPeakIndex(spectrum, fragments[f].Mz, tolerance, _parameters.MostIntenseFragmentPeak, bins.Hint(fragments[f].Mz));
                if (peakIndices is not null)
                    peakIndices[f][k] = i;
                if (i >= 0)
                {
                    traces[f][k] = spectrum.YArray[i];
                    ppm[f][k] = Math.Abs(spectrum.XArray[i] - fragments[f].Mz) / fragments[f].Mz * 1e6;
                }
            }
        }

        // The core: the most intense library fragments, in library rank; the rest enter as their own feature
        int[] coreIndices = FragmentCoElution.TopIndices(libraryIntensities, CoreFragmentCount)
            .OrderByDescending(i => libraryIntensities[i]).ThenBy(i => i).ToArray();
        var core = coreIndices.Select(i => traces[i]).ToList();
        var rest = new List<double[]>(fragments.Count);
        for (int i = 0; i < fragments.Count; i++)
            if (Array.IndexOf(coreIndices, i) < 0)
                rest.Add(traces[i]);
        // The core again, counting only peaks within a tight fraction of the tolerance
        double tightPpm = TightToleranceFraction * _parameters.FragmentTolerancePpm;
        var tightCore = new List<double[]>(coreIndices.Length);
        foreach (int i in coreIndices)
        {
            var tight = new double[reachable.Length];
            for (int k = 0; k < tight.Length; k++)
                tight[k] = ppm[i][k] <= tightPpm ? traces[i][k] : 0;
            tightCore.Add(tight);
        }
        // And within 0.2 x the tolerance, for each fragment's best-tolerance co-elution (MaxToleranceCoElution, DIA-NN's)
        List<double[]>? finestCore = null;
        if (_parameters.MaxToleranceCoElution)
        {
            double finestPpm = 0.2 * _parameters.FragmentTolerancePpm;
            finestCore = new List<double[]>(coreIndices.Length);
            foreach (int i in coreIndices)
            {
                var finest = new double[reachable.Length];
                for (int k = 0; k < finest.Length; k++)
                    finest[k] = ppm[i][k] <= finestPpm ? traces[i][k] : 0;
                finestCore.Add(finest);
            }
        }
        double[] summed = new double[reachable.Length];
        for (int k = 0; k < summed.Length; k++)
        {
            double sum = 0;
            foreach (var trace in traces)
                sum += trace[k];
            summed[k] = sum;
        }
        // Every scan's apex score, to judge how far a candidate stands out from the rest of its window (PECAN, OpenSWATH)
        double[] apexScores = FragmentCoElution.ApexScores(traces, libraryIntensities, _parameters.ApexHalfWidthScans);
        double[] scored = apexScores.Where(v => v > 0).ToArray();
        double scoreMean = scored.Length > 0 ? scored.Average() : 0;
        double scoreSd = scored.Length > 1 ? Math.Sqrt(scored.Sum(v => (v - scoreMean) * (v - scoreMean)) / (scored.Length - 1)) : 0;
        double windowSignal = summed.Sum();

        // The precursor's MS1 traces, monoisotopic and M+1, in the MS1 scan nearest each reachable MS2 scan
        double[]? ms1Mono = null, ms1Isotope = null, ms1MonoTight = null, ms1MonoPpm = null, ms1Isotope2 = null;
        if (_ms1Index is not null)
        {
            ms1Mono = new double[reachable.Length];
            ms1Isotope = new double[reachable.Length];
            // Ms1TightCorrelation: the M0 trace again, keeping only peaks within half the MS1 tolerance (the closest peak is
            // the one found, so any peak that close would be it)
            if (_parameters.Ms1TightCorrelation)
                ms1MonoTight = new double[reachable.Length];
            // Ms1PeakFeatures: each scan's M0 mass error, and the M+2 trace
            if (_parameters.Ms1PeakFeatures)
            {
                ms1MonoPpm = new double[reachable.Length];
                ms1Isotope2 = new double[reachable.Length];
            }
            double isotopeMz = candidate.PrecursorMz + 1.0033548 / Math.Max(1, (int)candidate.Charge);
            for (int k = 0; k < reachable.Length; k++)
            {
                double scanRt = scans[reachable[k]].RetentionTime;
                int ms1 = NearestMs1(scanRt);
                double monoMz = Ms1Target(candidate.PrecursorMz, scanRt);
                var mono = _ms1Index.GetIndexedPeak(monoMz, ms1, Ms1Tolerance);
                ms1Mono[k] = mono?.Intensity ?? 0;
                if (ms1MonoPpm is not null && ms1Isotope2 is not null)
                {
                    ms1MonoPpm[k] = mono is null ? 0 : Math.Abs(mono.M - monoMz) / monoMz * 1e6;
                    ms1Isotope2[k] = _ms1Index.GetIndexedPeak(Ms1Target(candidate.PrecursorMz + 2 * 1.0033548 / Math.Max(1, (int)candidate.Charge), scanRt), ms1, Ms1Tolerance)?.Intensity ?? 0;
                }
                if (ms1MonoTight is not null && mono is not null
                    && Math.Abs(mono.M - monoMz) / monoMz * 1e6 <= Ms1Tolerance.Value / 2)
                    ms1MonoTight[k] = mono.Intensity;
                ms1Isotope[k] = _ms1Index.GetIndexedPeak(Ms1Target(isotopeMz, scanRt), ms1, Ms1Tolerance)?.Intensity ?? 0;
            }
        }
        // The gate (MinimumApexFragments) picks candidates only; the apex-score features still see every scan
        double[] candidateScores = apexScores;
        if (_parameters.MinimumApexFragments > 0)
        {
            candidateScores = (double[])apexScores.Clone();
            for (int s = 0; s < candidateScores.Length; s++)
                if (core.Count(trace => trace[s] > 0) < _parameters.MinimumApexFragments)
                    candidateScores[s] = 0;
        }
        var apexes = _parameters.DiaNnPeakFinding
            ? DiaNnCandidateApexes(core, _parameters.ApexHalfWidthScans, _parameters.MaxApexCandidates, _parameters.DiaNnPeakFindingMs1 ? ms1Mono : null, _parameters.CandidatePeaksFromSummedCore).ToList()
            : FragmentCoElution.FindApexes(candidateScores, _parameters.ApexHalfWidthScans, _parameters.MaxApexCandidates).ToList();
        if (_parameters.FragmentApexCandidate && reachable.Length > 1)
        {
            // The best single fragment's own apex over the whole window, if no candidate is already there
            double[] bestTrace = FragmentCoElution.Smooth(core[FragmentCoElution.BestFragment(core, 0, reachable.Length - 1)]);
            int fragmentApex = 0;
            for (int s = 1; s < bestTrace.Length; s++)
                if (bestTrace[s] > bestTrace[fragmentApex])
                    fragmentApex = s;
            if (bestTrace[fragmentApex] > 0 && apexes.All(a => Math.Abs(a - fragmentApex) > _parameters.ApexHalfWidthScans))
                apexes.Add(fragmentApex);
        }
        // The extra fragments, the unfragmented precursor and the isotope shadows enter only through correlations over a
        // candidate's co-elution window (and its apex), so they are read only there: the same scores, far fewer lookups
        int coElutionHalfWidth = _parameters.CoElutionHalfWidthScans ?? _parameters.ApexHalfWidthScans;
        var inWindow = new bool[reachable.Length];
        foreach (int apex in apexes)
            for (int k = Math.Max(0, apex - coElutionHalfWidth); k <= Math.Min(reachable.Length - 1, apex + coElutionHalfWidth); k++)
                inWindow[k] = true;
        var extraTraces = extraFragments.Select(_ => new double[reachable.Length]).ToArray();
        if (extraTraces.Length > 0)
            for (int k = 0; k < reachable.Length; k++)
            {
                if (!inWindow[k])
                    continue;
                var spectrum = scans[reachable[k]].MassSpectrum;
            var bins = _binsByScan.GetOrAdd(scans[reachable[k]], s => PeakBins.Of(s.MassSpectrum));
                for (int x = 0; x < extraTraces.Length; x++)
                {
                    int i = FragmentPeakIndex(spectrum, extraFragments[x].Mz, tolerance, _parameters.MostIntenseFragmentPeak, bins.Hint(extraFragments[x].Mz));
                    if (i >= 0)
                        extraTraces[x][k] = spectrum.YArray[i];
                }
            }

        // DIA-NN's extra traces: the unfragmented precursor m/z in MS2, and one isotope below each core fragment
        double[]? unfragmented = null;
        List<double[]>? shadows = null, spikeSuppressed = null;
        if (_parameters.DiaNnScores)
        {
            unfragmented = new double[reachable.Length];
            shadows = coreIndices.Select(_ => new double[reachable.Length]).ToList();
            spikeSuppressed = core.Select(SpikeSuppressed).ToList();
            for (int k = 0; k < reachable.Length; k++)
            {
                if (!inWindow[k])
                    continue;
                var spectrum = scans[reachable[k]].MassSpectrum;
            var bins = _binsByScan.GetOrAdd(scans[reachable[k]], s => PeakBins.Of(s.MassSpectrum));
                int i = FragmentPeakIndex(spectrum, candidate.PrecursorMz, tolerance, _parameters.MostIntenseFragmentPeak, bins.Hint(candidate.PrecursorMz));
                if (i >= 0)
                    unfragmented[k] = spectrum.YArray[i];
                for (int c = 0; c < coreIndices.Length; c++)
                {
                    var fragment = fragments[coreIndices[c]];
                    double shadowMz = fragment.Mz - 1.0033548 / Math.Max(1, fragment.Charge);
                    int found = peakIndices![coreIndices[c]][k];
                    i = found >= 0 ? FragmentPeakIndex(spectrum, shadowMz, tolerance, _parameters.MostIntenseFragmentPeak, found)
                        : FragmentPeakIndex(spectrum, shadowMz, tolerance, _parameters.MostIntenseFragmentPeak, bins.Hint(shadowMz));
                    if (i >= 0)
                        shadows[c][k] = spectrum.YArray[i];
                }
            }
        }

        foreach (int apex in apexes)
        {
            double[] apexIntensities = traces.Select(trace => trace[apex]).ToArray();
            // Co-elution against the best of the six most intense library fragments, smoothed: one reliable profile rather
            // than an average that an interfered fragment drags along (DIA-NN's approach, Demichev et al. 2020)
            int from = Math.Max(0, apex - coElutionHalfWidth);
            int to = Math.Min(reachable.Length - 1, apex + coElutionHalfWidth);
            double coElution = 0, tightCoElution = 0, remainingCoElution = 0, maxToleranceCoElution = 0, sqrtCoElution = 0;
            var fragmentCorrelations = new double[CoreFragmentCount];
            double[]? reference = null;
            if (to > from)
            {
                reference = FragmentCoElution.Smooth(core[FragmentCoElution.BestFragment(core, from, to)]);
                double[] correlations = FragmentCoElution.CorrelationsTo(core, reference, from, to);
                coElution = correlations.Average();
                Array.Copy(correlations, fragmentCorrelations, correlations.Length);
                if (rest.Count > 0)
                    remainingCoElution = FragmentCoElution.CorrelationsTo(rest, reference, from, to).Average();
                var tightReference = FragmentCoElution.Smooth(tightCore[FragmentCoElution.BestFragment(tightCore, from, to)]);
                tightCoElution = FragmentCoElution.CorrelationsTo(tightCore, tightReference, from, to).Average();
                if (finestCore is not null)
                {
                    double[] atTight = FragmentCoElution.CorrelationsTo(tightCore, reference, from, to);
                    double[] atFinest = FragmentCoElution.CorrelationsTo(finestCore, reference, from, to);
                    double sum = 0;
                    for (int f = 0; f < correlations.Length; f++)
                        sum += Math.Max(correlations[f], Math.Max(atTight[f], atFinest[f]));
                    maxToleranceCoElution = sum / correlations.Length;
                }
                // The same on square-root traces: the apex scans weigh less against the flanks
                if (_parameters.SqrtCoElution)
                {
                    double[][] sqrtCore = core.Select(trace => trace.Select(Math.Sqrt).ToArray()).ToArray();
                    sqrtCoElution = FragmentCoElution.CorrelationsTo(sqrtCore, reference.Select(Math.Sqrt).ToArray(), from, to).Average();
                }
            }
            double cosine = SpectralSimilarity.CosineOfAlignedVectors(apexIntensities, libraryIntensities);

            // Library similarity across the peak, each scan weighted by the profile squared
            double windowCosine = cosine;
            var observedAt = new double[traces.Length];
            if (reference is not null)
            {
                double weighted = 0, weights = 0;
                for (int s = from; s <= to; s++)
                {
                    double weight = reference[s] * reference[s];
                    if (weight <= 0)
                        continue;
                    for (int f = 0; f < traces.Length; f++)
                        observedAt[f] = traces[f][s];
                    weighted += weight * SpectralSimilarity.CosineOfAlignedVectors(observedAt, libraryIntensities);
                    weights += weight;
                }
                if (weights > 0)
                    windowCosine = weighted / weights;
            }

            // Mass accuracy of the fragments seen at the apex, from the peaks their intensities were read from
            var apexSpectrum = scans[reachable[apex]].MassSpectrum;
            var ppmErrors = new List<double>();
            for (int f = 0; f < fragments.Count; f++)
                if (apexIntensities[f] > 0)
                    ppmErrors.Add(ppm[f][apex]);

            var (peakStart, peakEnd) = FragmentCoElution.PeakBounds(summed, apex);

            // Median-profile co-elution (EncyclopeDIA, AlphaDIA): each core trace normalised to sum 1 over the window,
            // the per-scan median as a reference that no single fragment can hijack
            double medianCoElution = 0, goodFragments = 0, veryGoodFragment = 0, gaussianFit = 0;
            if (to > from)
            {
                double[] median = MedianProfile(core, from, to, reachable.Length);
                double[] toMedian = FragmentCoElution.CorrelationsTo(core, median, from, to);
                medianCoElution = toMedian.Average();
                goodFragments = toMedian.Count(r => r >= 0.75);
                veryGoodFragment = toMedian.Any(r => r >= 0.9) ? 1 : 0;
                double sigma = Math.Max(1, (to - from) / 4.0);
                double[] gaussian = Enumerable.Range(0, reachable.Length).Select(s => Math.Exp(-0.5 * Math.Pow((s - apex) / sigma, 2))).ToArray();
                gaussianFit = FragmentCoElution.CorrelationsTo([median], gaussian, from, to)[0];
            }

            // Library agreement on peak areas rather than one apex scan (OpenSWATH, Skyline): cosine of square-root areas,
            // Pearson of areas, and Manhattan distance of sum-normalised areas
            double[] areas = traces.Select(trace => Area(trace, peakStart, peakEnd)).ToArray();
            double areaSqrtCosine = SpectralSimilarity.CosineOfAlignedVectors(areas.Select(Math.Sqrt).ToArray(), libraryIntensities.Select(Math.Sqrt).ToArray());
            double areaPearson = areas.Any(a => a > 0) ? Math.Max(0, FragmentCoElution.CorrelationsTo([areas], libraryIntensities, 0, areas.Length - 1)[0]) : 0;
            double areaSum = areas.Sum(), librarySum = libraryIntensities.Sum();
            double areaManhattan = areaSum > 0 ? areas.Select((a, f) => Math.Abs(a / areaSum - libraryIntensities[f] / librarySum)).Sum() : 2;

            // Fewer fragments, weighted by expectation: a weak precursor shows only its strongest fragments above noise, so
            // score it on those rather than penalise it for fragments it could not show
            double[] allCorrelations = reference is null ? new double[fragments.Count] : FragmentCoElution.CorrelationsTo(traces, reference, from, to);
            double libraryTotal = libraryIntensities.Sum();
            double weightedCoElution = libraryIntensities.Select((w, f) => w * allCorrelations[f]).Sum() / libraryTotal;
            double weightedMatchedFraction = libraryIntensities.Select((w, f) => apexIntensities[f] > 0 ? w : 0).Sum() / libraryTotal;
            int[] byRank = Enumerable.Range(0, fragments.Count).OrderByDescending(f => libraryIntensities[f]).ThenBy(f => f).ToArray();
            int[] top3 = byRank.Take(3).ToArray();
            double top3CoElution = top3.Average(f => allCorrelations[f]);
            double top3Cosine = SpectralSimilarity.CosineOfAlignedVectors(top3.Select(f => apexIntensities[f]).ToArray(), top3.Select(f => libraryIntensities[f]).ToArray());
            double top1Present = apexIntensities[byRank[0]] > 0 ? 1 : 0;
            double top2Present = byRank.Length > 1 && apexIntensities[byRank[1]] > 0 ? 1 : 0;
            double top1Ppm = top1Present > 0 ? ppm[byRank[0]][apex] : _parameters.FragmentTolerancePpm;

            // Detectable fragments: the library pattern scaled to the apex by least squares, against the apex scan's noise
            // (its 25th-percentile intensity); fragments expected below 3x noise are not held against the precursor
            double libraryScale = libraryIntensities.Select((l, f) => l * apexIntensities[f]).Sum() / libraryIntensities.Sum(l => l * l);
            double noise = _noiseByScan.GetOrAdd(scans[reachable[apex]], NoiseOf);
            int[] detectable = Enumerable.Range(0, fragments.Count).Where(f => libraryScale * libraryIntensities[f] >= 3 * noise).ToArray();
            if (detectable.Length == 0)
                detectable = [byRank[0]];
            double detectableMatchedFraction = detectable.Count(f => apexIntensities[f] > 0) / (double)detectable.Length;
            double detectableCoElution = detectable.Average(f => allCorrelations[f]);

            // Fragment isotopes at the apex (OpenSWATH): a matched peak with a bigger peak one isotope below is probably another
            // ion's isotope; a real fragment usually shows its own M+1
            double isotopeOverlap = 0, fragmentM1 = 0, matchedWeight = 0;
            for (int f = 0; f < fragments.Count; f++)
            {
                if (apexIntensities[f] <= 0)
                    continue;
                double w = libraryIntensities[f];
                matchedWeight += w;
                double spacing = 1.0033548 / Math.Max(1, fragments[f].Charge);
                int below = apexSpectrum.GetClosestPeakIndex(fragments[f].Mz - spacing);
                if (tolerance.Within(apexSpectrum.XArray[below], fragments[f].Mz - spacing) && apexSpectrum.YArray[below] > apexIntensities[f])
                    isotopeOverlap += w;
                int above = apexSpectrum.GetClosestPeakIndex(fragments[f].Mz + spacing);
                if (tolerance.Within(apexSpectrum.XArray[above], fragments[f].Mz + spacing))
                    fragmentM1 += w;
            }
            if (matchedWeight > 0)
            {
                isotopeOverlap /= matchedWeight;
                fragmentM1 /= matchedWeight;
            }

            // Specificity (PECAN's idea): a match where peaks are sparse is less likely by chance than one in the crowded
            // low-m/z region. Chance probability p = 1 - exp(-density x window), with density the apex spectrum's peaks per
            // Th within +-5 Th and window the tolerance width in Th; each fragment's evidence is its surprisal -log10 p
            double matchSurprisal = 0, possibleSurprisal = 0, specificityWeightedCoElution = 0;
            for (int f = 0; f < fragments.Count; f++)
            {
                double mz = fragments[f].Mz;
                int lo = apexSpectrum.Size == 0 ? 0 : LowerBound(apexSpectrum.XArray, mz - 5);
                int hi = apexSpectrum.Size == 0 ? 0 : LowerBound(apexSpectrum.XArray, mz + 5);
                double density = Math.Max(1, hi - lo) / 10.0;
                double window = 2 * mz * _parameters.FragmentTolerancePpm * 1e-6;
                double chance = Math.Clamp(1 - Math.Exp(-density * window), 1e-6, 1);
                double surprisal = -Math.Log10(chance);
                possibleSurprisal += surprisal;
                specificityWeightedCoElution += surprisal * allCorrelations[f];
                if (apexIntensities[f] > 0)
                    matchSurprisal += surprisal;
            }
            double surprisalFraction = possibleSurprisal > 0 ? matchSurprisal / possibleSurprisal : 0;
            specificityWeightedCoElution = possibleSurprisal > 0 ? specificityWeightedCoElution / possibleSurprisal : 0;

            // MS1 as DIA-NN uses it: the precursor's MS1 trace (monoisotopic and M+1) correlated with the fragment profile,
            // never its raw intensity
            double ms1Correlation = 0, ms1IsotopeCorrelation = 0, ms1TightCorrelation = 0, ms1PeakMassError = 0, ms1Isotope2Correlation = 0;
            if (reference is not null && ms1Mono is not null && ms1Isotope is not null)
            {
                int ms1From = from, ms1To = to;
                if (_parameters.Ms1CorrelationHalfWidthScans is int ms1HalfWidth)
                {
                    ms1From = Math.Max(0, apex - ms1HalfWidth);
                    ms1To = Math.Min(reachable.Length - 1, apex + ms1HalfWidth);
                }
                ms1Correlation = FragmentCoElution.CorrelationsTo([ms1Mono], reference, ms1From, ms1To)[0];
                if (ms1MonoTight is not null)
                    ms1TightCorrelation = FragmentCoElution.CorrelationsTo([ms1MonoTight], reference, from, to)[0];
                if (ms1MonoPpm is not null && ms1Isotope2 is not null)
                {
                    // Intensity-weighted M0 mass error over the window; the full tolerance when no M0 was found there
                    double weighted = 0, weights = 0;
                    for (int s = from; s <= to; s++)
                        if (ms1Mono[s] > 0)
                        {
                            weighted += ms1Mono[s] * ms1MonoPpm[s];
                            weights += ms1Mono[s];
                        }
                    ms1PeakMassError = weights > 0 ? weighted / weights : Ms1Tolerance.Value;
                    ms1Isotope2Correlation = FragmentCoElution.CorrelationsTo([ms1Isotope2], reference, from, to)[0];
                }
                ms1IsotopeCorrelation = FragmentCoElution.CorrelationsTo([ms1Isotope], reference, ms1From, ms1To)[0];
            }

            // MS1 evidence that stands without the fragments, at the MS1 scan nearest the apex: the M0-M3 envelope against the
            // expected isotope pattern, the M0 mass error, and the apex's share of the window's MS1 trace maximum. On the
            // whole-proteome library a fifth of DIA-NN's identifications we missed had a clean MS1 peak and faint fragments.
            double ms1EnvelopeCosine = 0, ms1MassErrorPpm = Ms1Tolerance.Value, ms1ApexShare = 0, ms1ApexSignedErrorPpm = double.NaN;
            if (_ms1Index is not null && ms1Mono is not null)
            {
                double apexScanRt = scans[reachable[apex]].RetentionTime;
                int ms1Scan = NearestMs1(apexScanRt);
                int z = Math.Max(1, (int)candidate.Charge);
                double[] expected = ExpectedIsotopes((candidate.PrecursorMz - 1.007276) * z, 4);
                var observed = new double[4];
                for (int k = 0; k < 4; k++)
                {
                    double target = Ms1Target(candidate.PrecursorMz + k * 1.0033548 / z, apexScanRt);
                    var peak = _ms1Index.GetIndexedPeak(target, ms1Scan, Ms1Tolerance);
                    observed[k] = peak?.Intensity ?? 0;
                    if (k == 0 && peak is not null)
                    {
                        ms1MassErrorPpm = Math.Abs(peak.M - target) / target * 1e6;
                        ms1ApexSignedErrorPpm = (peak.M - candidate.PrecursorMz) / candidate.PrecursorMz * 1e6;
                    }
                }
                ms1EnvelopeCosine = observed.Any(v => v > 0) ? SpectralSimilarity.CosineOfAlignedVectors(observed, expected) : 0;
                double windowMax = ms1Mono.Max();
                double atApex = ms1Mono.Skip(Math.Max(0, apex - 1)).Take(3).Max();
                ms1ApexShare = windowMax > 0 ? atApex / windowMax : 0;
            }

            // Ms1PeakEnvelope: M-1 to M3 summed over every MS1 scan of the co-elution window, so a faint precursor's envelope
            // rests on several scans rather than the one nearest the apex
            // Ms1PeakEnvelopeTight: the same from peaks within 0.6x the MS1 tolerance, so the envelope is seen at two widths
            double ms1PeakEnvelopeCosine = 0, ms1PeakEnvelopeTightCosine = 0;
            if ((_parameters.Ms1PeakEnvelope || _parameters.Ms1PeakEnvelopeTight) && _ms1Index is not null)
            {
                int z = Math.Max(1, (int)candidate.Charge);
                var summedEnvelope = new double[5];
                var tightEnvelope = new double[5];
                double envelopeTightPpm = 0.6 * Ms1Tolerance.Value;
                int lastMs1 = -1;
                for (int s = from; s <= to; s++)
                {
                    double scanRt = scans[reachable[s]].RetentionTime;
                    int ms1Scan = NearestMs1(scanRt);
                    if (ms1Scan == lastMs1)
                        continue;
                    lastMs1 = ms1Scan;
                    for (int k = -1; k <= 3; k++)
                    {
                        double target = Ms1Target(candidate.PrecursorMz + k * 1.0033548 / z, scanRt);
                        var peak = _ms1Index.GetIndexedPeak(target, ms1Scan, Ms1Tolerance);
                        if (peak is null)
                            continue;
                        summedEnvelope[k + 1] += peak.Intensity;
                        if (Math.Abs(peak.M - target) / target * 1e6 <= envelopeTightPpm)
                            tightEnvelope[k + 1] += peak.Intensity;
                    }
                }
                double[] expectedEnvelope = ExpectedIsotopes((candidate.PrecursorMz - 1.007276) * z, 4);
                if (_parameters.Ms1PeakEnvelope)
                    ms1PeakEnvelopeCosine = PeakEnvelopeCosine(summedEnvelope, expectedEnvelope);
                if (_parameters.Ms1PeakEnvelopeTight)
                    ms1PeakEnvelopeTightCosine = PeakEnvelopeCosine(tightEnvelope, expectedEnvelope);
            }

            // Ms1Intensity: how much precursor there is over the peak (each MS1 scan once), and against the fragments' area
            double ms1LogIntensity = 0, ms1ToFragmentLogRatio = 0;
            if (_parameters.Ms1Intensity && ms1Mono is not null)
            {
                double ms1Sum = 0;
                int lastMs1 = -1;
                for (int s = peakStart; s <= peakEnd; s++)
                {
                    int ms1Scan = NearestMs1(scans[reachable[s]].RetentionTime);
                    if (ms1Scan == lastMs1)
                        continue;
                    lastMs1 = ms1Scan;
                    ms1Sum += ms1Mono[s];
                }
                ms1LogIntensity = Math.Log10(1 + ms1Sum);
                ms1ToFragmentLogRatio = ms1LogIntensity - Math.Log10(1 + areaSum);
            }

            // Ms1PeakSignalFraction: the share of the window's M0 signal inside this peak, as PeakSignalFraction does for fragments
            double ms1PeakSignalFraction = 0;
            if (_parameters.Ms1PeakSignalFraction && ms1Mono is not null)
            {
                double ms1Window = ms1Mono.Sum();
                double ms1Peak = 0;
                for (int s = peakStart; s <= peakEnd; s++)
                    ms1Peak += ms1Mono[s];
                ms1PeakSignalFraction = ms1Window > 0 ? ms1Peak / ms1Window : 0;
            }

            // Ms1PeakPoints: in how many of the peak's distinct MS1 scans the precursor is seen at all
            double ms1PeakPointFraction = 0;
            if (_parameters.Ms1PeakPoints && ms1Mono is not null)
            {
                int ms1Scans = 0, ms1Seen = 0, lastMs1 = -1;
                for (int s = peakStart; s <= peakEnd; s++)
                {
                    int ms1Scan = NearestMs1(scans[reachable[s]].RetentionTime);
                    if (ms1Scan == lastMs1)
                        continue;
                    lastMs1 = ms1Scan;
                    ms1Scans++;
                    if (ms1Mono[s] > 0)
                        ms1Seen++;
                }
                ms1PeakPointFraction = ms1Scans > 0 ? (double)ms1Seen / ms1Scans : 0;
            }

            double extraCoElution = 0, extraMatchedFraction = 0, extraWeightedCoElution = 0;
            if (reference is not null && extraTraces.Length > 0)
            {
                double[] correlations = FragmentCoElution.CorrelationsTo(extraTraces, reference, from, to);
                extraCoElution = correlations.Average();
                extraMatchedFraction = extraTraces.Count(trace => trace[apex] > 0) / (double)extraTraces.Length;
                double weight = extraFragments.Sum(LibraryIntensity);
                extraWeightedCoElution = weight > 0 ? extraFragments.Select((f, x) => LibraryIntensity(f) * correlations[x]).Sum() / weight : 0;
            }
            else if (_parameters.ExtraFragmentFallback && extraTraces.Length == 0)
            {
                // No fragments beyond the top N: repeat the main fragments' evidence rather than zeros that read as failure
                extraCoElution = coElution;
                extraMatchedFraction = (double)ppmErrors.Count / fragments.Count;
                extraWeightedCoElution = weightedCoElution;
            }

            double minCorr = 0, nfCorr = 0, shadowCorr = 0;
            if (reference is not null && unfragmented is not null && shadows is not null && spikeSuppressed is not null)
            {
                minCorr = FragmentCoElution.CorrelationsTo(spikeSuppressed, reference, from, to).Average();
                nfCorr = FragmentCoElution.CorrelationsTo([unfragmented], reference, from, to)[0];
                shadowCorr = FragmentCoElution.CorrelationsTo(shadows, reference, from, to).Average();
            }

            // DIA-NN's pSig: each core fragment's share of the core's summed signal across the co-elution window, in library rank
            var signalShares = new double[CoreFragmentCount];
            if (_parameters.DiaNnSignalShare)
            {
                double total = 0;
                for (int c = 0; c < core.Count; c++)
                {
                    for (int s = from; s <= to; s++)
                        signalShares[c] += core[c][s];
                    total += signalShares[c];
                }
                for (int c = 0; c < signalShares.Length; c++)
                    signalShares[c] = total > 0 ? signalShares[c] / total : 0;
            }
            // Uniqueness: this apex against the best competing scan outside its co-elution window, its z-score among the
            // window's scans, and the share of the window's fragment signal inside its peak
            double apexScore = apexScores[apex];
            double competitor = 0;
            bool anyCompetitor = false;
            for (int s = 0; s < apexScores.Length; s++)
                if (Math.Abs(s - apex) > _parameters.ApexHalfWidthScans)
                {
                    competitor = anyCompetitor ? Math.Max(competitor, apexScores[s]) : apexScores[s];
                    anyCompetitor = true;
                }
            double apexScoreDelta = apexScore > 0 ? Math.Clamp((apexScore - competitor) / apexScore, -1, 1) : 0;
            double apexScoreZ = scoreSd > 0 ? (apexScore - scoreMean) / scoreSd : 0;
            double peakSignalFraction = windowSignal > 0 ? Enumerable.Range(peakStart, peakEnd - peakStart + 1).Sum(s => summed[s]) / windowSignal : 0;
            var apexRt = new RtMinutes(scans[reachable[apex]].RetentionTime);
            var apexIrt = _irtMap.ToIrt(apexRt);

            double[] features =
            [
                cosine,
                coElution,
                (double)ppmErrors.Count / fragments.Count,
                ppmErrors.Count > 0 ? ppmErrors.Average() : _parameters.FragmentTolerancePpm,
                Math.Abs(apexIrt.Value - expectedIrt),
                Math.Log10(1 + apexIntensities.Sum()),
                peakEnd - peakStart + 1,
                .. fragmentCorrelations,
                tightCoElution,
                remainingCoElution,
                windowCosine,
                medianCoElution,
                goodFragments,
                veryGoodFragment,
                gaussianFit,
                areaSqrtCosine,
                areaPearson,
                areaManhattan,
                apexScoreDelta,
                apexScoreZ,
                peakSignalFraction,
                weightedCoElution,
                weightedMatchedFraction,
                top3CoElution,
                top3Cosine,
                top1Present,
                top2Present,
                top1Ppm,
                detectableMatchedFraction,
                detectableCoElution,
                detectable.Length,
                ms1Correlation,
                ms1IsotopeCorrelation,
                isotopeOverlap,
                fragmentM1,
                matchSurprisal,
                surprisalFraction,
                specificityWeightedCoElution,
                entry.BaseSequence.Length,
                ms1EnvelopeCosine,
                ms1MassErrorPpm,
                ms1ApexShare,
                0,   // SiblingCoElution, filled in once every candidate is scored (AddSiblingFeatures)
                1.0, // SiblingApexDeltaMinutes, likewise
                candidate.PrecursorMz,
                (double)candidate.Charge,
                entry.MatchedFragmentIons.Count,
                minCorr,
                nfCorr,
                shadowCorr,
                .. signalShares,
                extraCoElution,
                extraMatchedFraction,
                extraWeightedCoElution,
                maxToleranceCoElution,
                ms1TightCorrelation,
                ms1PeakMassError,
                ms1Isotope2Correlation,
                sqrtCoElution,
                ms1PeakEnvelopeCosine,
                ms1PeakEnvelopeTightCosine,
                ms1LogIntensity,
                ms1ToFragmentLogRatio,
                ms1PeakSignalFraction,
                ms1PeakPointFraction,
                _parameters.LibraryIrtFeature ? candidate.Irt : 0,
            ];

            yield return new DiaPrecursorMatch(
                candidate.PrecursorIdx,
                entry.FullSequence,
                candidate.Charge,
                candidate.PrecursorMz,
                candidate.IsDecoy != 0,
                new Irt(candidate.Irt),
                apexRt,
                apexIrt,
                cosine * coElution,
                features,
                Quantity: areaSum) { Ms1ApexErrorPpm = ms1ApexSignedErrorPpm, ApexScanNumber = scans[reachable[apex]].OneBasedScanNumber };
        }
    }

    /// <summary>
    /// The peak a fragment is read from: the one nearest <paramref name="mz"/>, or with <paramref name="mostIntense"/> the most
    /// intense within tolerance (DIA-NN 1.8's level()); -1 when none is within tolerance.
    /// </summary>
    public static int FragmentPeakIndex(MzSpectrum spectrum, double mz, Tolerance tolerance, bool mostIntense) =>
        spectrum.Size == 0 ? -1 : PeakFromNearest(spectrum, mz, tolerance, mostIntense, spectrum.GetClosestPeakIndex(mz));

    /// <summary>
    /// <see cref="FragmentPeakIndex(MzSpectrum, double, Tolerance, bool)"/> found by walking from <paramref name="hint"/>, an
    /// index near <paramref name="mz"/>, rather than by binary search: the same peak, cheaper when the hint is a few peaks
    /// away (a fragment's isotope shadow, one isotope below the fragment just read).
    /// </summary>
    public static int FragmentPeakIndex(MzSpectrum spectrum, double mz, Tolerance tolerance, bool mostIntense, int hint)
    {
        int n = spectrum.Size;
        if (n == 0)
            return -1;
        var x = spectrum.XArray;
        // The first index at or above mz, then GetClosestIndex's rule: an exact match, an end, or the nearer neighbour
        // (the upper one on a tie)
        int i = Math.Clamp(hint, 0, n - 1);
        while (i < n && x[i] < mz) i++;
        while (i > 0 && x[i - 1] >= mz) i--;
        int nearest = i == n ? n - 1 : i == 0 || x[i] == mz ? i : mz - x[i - 1] < x[i] - mz ? i - 1 : i;
        return PeakFromNearest(spectrum, mz, tolerance, mostIntense, nearest);
    }

    /// <summary>
    /// A spectrum's m/z bins (<see cref="Width"/> Th): the first peak at or above each bin's lower edge, the start from which
    /// <see cref="FragmentPeakIndex(MzSpectrum, double, Tolerance, bool, int)"/> walks to a fragment's peak instead of a binary
    /// search (the largest share of extraction, with 25 lookups per scan per precursor). Built once per scan.
    /// </summary>
    public sealed class PeakBins
    {
        public const double Width = 2;
        private readonly int _lowBin;
        private readonly int[] _starts;
        private readonly int _last;

        private PeakBins(int lowBin, int[] starts, int last)
        {
            _lowBin = lowBin;
            _starts = starts;
            _last = last;
        }

        public static PeakBins Of(MzSpectrum spectrum)
        {
            var x = spectrum.XArray;
            if (x.Length == 0)
                return new PeakBins(0, [0], 0);
            int low = (int)Math.Floor(x[0] / Width), high = (int)Math.Floor(x[^1] / Width);
            var starts = new int[high - low + 1];
            int i = 0;
            for (int b = 0; b < starts.Length; b++)
            {
                double edge = (low + b) * Width;
                while (i < x.Length && x[i] < edge) i++;
                starts[b] = Math.Min(i, x.Length - 1);
            }
            return new PeakBins(low, starts, x.Length - 1);
        }

        /// <summary>A valid index near <paramref name="mz"/>: the first peak of its bin, clamped at either end.</summary>
        public int Hint(double mz)
        {
            double bin = Math.Floor(mz / Width) - _lowBin;
            return bin < 0 ? 0 : bin >= _starts.Length ? _last : _starts[(int)bin];
        }
    }

    private readonly System.Collections.Concurrent.ConcurrentDictionary<MsDataScan, PeakBins> _binsByScan = new(ReferenceEqualityComparer.Instance);

    /// <summary>
    /// DIA-NN's fragment rule (MinFrAAs = 3, MinFrMz = 200, MaxFrMz = 1800): a terminal ion of at least 3 residues within
    /// 200-1800 m/z. Internal and diagnostic ions are never scored: FragmentNumber is an internal ion's start residue and a
    /// diagnostic ion has no residues (MSL thread, 010).
    /// </summary>
    public static bool IsDiaNnScorable(MslFragmentIon ion, double minimumMz = 200, int minimumResidues = 3, int maximumCharge = 0,
        int precursorCharge = 0, bool chargeBelowPrecursor = false) =>
        !ion.IsInternalFragment && !ion.IsDiagnosticIon && ion.FragmentNumber >= minimumResidues && ion.Mz >= minimumMz && ion.Mz <= 1800 && (maximumCharge <= 0 || ion.Charge <= maximumCharge)
        && (!chargeBelowPrecursor || precursorCharge <= 1 || ion.Charge < precursorCharge);

    private static int PeakFromNearest(MzSpectrum spectrum, double mz, Tolerance tolerance, bool mostIntense, int nearest)
    {
        if (!tolerance.Within(spectrum.XArray[nearest], mz))
            return -1;
        if (!mostIntense)
            return nearest;
        int best = nearest;
        for (int j = nearest - 1; j >= 0 && tolerance.Within(spectrum.XArray[j], mz); j--)
            if (spectrum.YArray[j] > spectrum.YArray[best]) best = j;
        for (int j = nearest + 1; j < spectrum.Size && tolerance.Within(spectrum.XArray[j], mz); j++)
            if (spectrum.YArray[j] > spectrum.YArray[best]) best = j;
        return best;
    }

    /// <summary>
    /// Candidate apexes as DIA-NN 1.8 finds them (Searcher::peaks). The score of a scan is its best fragment's summed
    /// correlation to the other fragments in <paramref name="core"/>, over the <paramref name="halfWidth"/> scans on either
    /// side. A scan is a candidate when at least two fragments are present, the score is at least 0.5, and the 1-2-1 smoothed
    /// best fragment peaks there (within a third of the half width). Candidates within 2.0 of the best are kept, more than
    /// <paramref name="halfWidth"/> apart, best first, at most <paramref name="maxCount"/>. DIA-NN's MS1 correlation term is
    /// added when <paramref name="ms1"/> is given: the MS1 trace's correlation with the best fragment. With
    /// <paramref name="summedCorePeaks"/>, scans where the sum of the smoothed fragments peaks also compete on score for the
    /// <paramref name="maxCount"/> places, if more than <paramref name="halfWidth"/> from every candidate, so none of those moves.
    /// </summary>
    public static int[] DiaNnCandidateApexes(IReadOnlyList<double[]> core, int halfWidth, int maxCount, double[]? ms1 = null, bool summedCorePeaks = false)
    {
        int n = core.Count == 0 ? 0 : core[0].Length;
        var smoothed = core.Select(FragmentCoElution.Smooth).ToArray();
        double[]? summed = null;
        if (summedCorePeaks)
        {
            summed = new double[n];
            foreach (var trace in smoothed)
                for (int s = 0; s < n; s++)
                    summed[s] += trace[s];
        }
        int local = Math.Max(1, halfWidth / 3);
        var candidates = new List<(int Scan, double Score)>();
        var summedOnly = new List<(int Scan, double Score)>();
        // Flat rather than [,]: a multi-dimensional array per call was a measurable share of extraction (CreateInstanceMDArray)
        var correlation = new double[core.Count * core.Count];
        for (int s = 0; s < n; s++)
        {
            int present = 0;
            foreach (var trace in core)
                if (trace[s] > 0) present++;
            if (present < 2)
                continue;
            int from = Math.Max(0, s - halfWidth), to = Math.Min(n - 1, s + halfWidth);
            for (int f = 0; f < core.Count; f++)
                for (int g = f + 1; g < core.Count; g++)
                    correlation[f * core.Count + g] = correlation[g * core.Count + f] = ClippedPearson(core[f], core[g], from, to);
            int best = 0;
            double bestSum = double.NegativeInfinity;
            for (int f = 0; f < core.Count; f++)
            {
                double sum = 0;
                for (int g = 0; g < core.Count; g++)
                    if (g != f) sum += correlation[f * core.Count + g];
                if (sum > bestSum) { bestSum = sum; best = f; }
            }
            if (ms1 is not null)
                bestSum += ClippedPearson(core[best], ms1, from, to);
            if (bestSum < 0.5)
                continue;
            if (IsLocalMaximum(smoothed[best], s, local))
                candidates.Add((s, bestSum));
            else if (summed is not null && IsLocalMaximum(summed, s, local))
                summedOnly.Add((s, bestSum));
        }
        var chosen = new List<(int Scan, double Score)>();
        double top = candidates.Count > 0 ? candidates.Max(c => c.Score) : 0;
        foreach (var (scan, score) in candidates.OrderByDescending(c => c.Score).ThenBy(c => c.Scan))
        {
            // With sum-peaks still to weigh, the cap waits until they have competed for its places
            if (score < top - 2.0 || (summedOnly.Count == 0 && chosen.Count == maxCount))
                break;
            if (chosen.All(c => Math.Abs(c.Scan - scan) > halfWidth))
                chosen.Add((scan, score));
        }
        if (summedOnly.Count == 0)
            return chosen.Select(c => c.Scan).ToArray();

        // A sum-peak can never sit on a candidate of the best-fragment rule: taking over a nearby apex moved found precursors off their
        // peaks and cost more than it rescued. But it does compete on score for the capped places, since the peaks it rescues mostly
        // belong to precursors whose places are already full of weaker candidates.
        double bar = (candidates.Count > 0 ? top : summedOnly.Max(c => c.Score)) - 2.0;
        foreach (var (scan, score) in summedOnly.OrderByDescending(c => c.Score).ThenBy(c => c.Scan))
        {
            if (score < bar)
                break;
            if (chosen.All(c => Math.Abs(c.Scan - scan) > halfWidth))
                chosen.Add((scan, score));
        }
        return chosen.OrderByDescending(c => c.Score).ThenBy(c => c.Scan).Take(maxCount).Select(c => c.Scan).ToArray();
    }

    /// <summary>Whether no scan within <paramref name="radius"/> of <paramref name="s"/> is higher; on a tie, the earliest scan is the maximum.</summary>
    private static bool IsLocalMaximum(double[] trace, int s, int radius)
    {
        for (int t = Math.Max(0, s - radius); t <= Math.Min(trace.Length - 1, s + radius); t++)
            if (trace[t] > trace[s] || (t < s && trace[t] == trace[s]))
                return false;
        return true;
    }

    /// <summary>
    /// A trace with single-scan spikes suppressed, as DIA-NN's pMinCorr reads it: each scan the minimum of itself and its
    /// neighbours (an end scan has one neighbour).
    /// </summary>
    public static double[] SpikeSuppressed(double[] trace)
    {
        var result = new double[trace.Length];
        for (int s = 0; s < trace.Length; s++)
        {
            double min = trace[s];
            if (s > 0) min = Math.Min(min, trace[s - 1]);
            if (s < trace.Length - 1) min = Math.Min(min, trace[s + 1]);
            result[s] = min;
        }
        return result;
    }

    /// <summary>Pearson's r over [from, to]; negative, flat or undefined counts as 0.</summary>
    internal static double ClippedPearson(double[] a, double[] b, int from, int to)
    {
        int points = to - from + 1;
        if (points < 3)
            return 0;
        double meanA = 0, meanB = 0;
        for (int i = from; i <= to; i++) { meanA += a[i]; meanB += b[i]; }
        meanA /= points;
        meanB /= points;
        double cov = 0, varA = 0, varB = 0;
        for (int i = from; i <= to; i++)
        {
            double da = a[i] - meanA, db = b[i] - meanB;
            cov += da * db;
            varA += da * da;
            varB += db * db;
        }
        double r = varA > 0 && varB > 0 ? cov / Math.Sqrt(varA * varB) : 0;
        return double.IsFinite(r) && r > 0 ? r : 0;
    }

    /// <summary>The first index whose value is at least <paramref name="x"/> in an ascending array.</summary>
    private static int LowerBound(double[] sorted, double x)
    {
        int lo = 0, hi = sorted.Length;
        while (lo < hi)
        {
            int mid = (lo + hi) / 2;
            if (sorted[mid] < x) lo = mid + 1; else hi = mid;
        }
        return lo;
    }

    /// <summary>
    /// Charge-state sibling features (AlphaDIA's elution groups): for each candidate, the best co-elution among candidates of
    /// the same sequence and decoy flag at another charge whose apex lies within <paramref name="rtToleranceMinutes"/>, and
    /// the RT gap to the nearest of them (capped at 1 minute). Only the sequence and decoy flag group them, so targets and
    /// decoys are treated alike.
    /// </summary>
    private static void AddSiblingFeatures(List<DiaPrecursorMatch> rows, double rtToleranceMinutes, bool topCandidateOnly = false)
    {
        int coElution = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "CoElution");
        int siblingCoElution = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "SiblingCoElution");
        int siblingDelta = Array.IndexOf(DiaPrecursorMatch.FeatureNames, "SiblingApexDeltaMinutes");
        foreach (var group in rows.GroupBy(r => (r.FullSequence, r.IsDecoy)))
        {
            var members = group.ToList();
            if (members.Select(r => r.Charge).Distinct().Count() < 2)
                continue;
            // Support from each charge state's top candidate (by the pre-rescoring score), or from any of its candidates
            var supporters = topCandidateOnly
                ? members.GroupBy(r => r.Charge).Select(g => g.OrderByDescending(r => r.Score).ThenBy(r => r.ApexRt.Value).First()).ToList()
                : members;
            foreach (var row in members)
            {
                double best = 0, nearest = 1.0;
                foreach (var other in supporters)
                {
                    if (other.Charge == row.Charge)
                        continue;
                    double gap = Math.Abs(other.ApexRt.Value - row.ApexRt.Value);
                    nearest = Math.Min(nearest, gap);
                    if (gap <= rtToleranceMinutes)
                        best = Math.Max(best, other.Features[coElution]);
                }
                row.Features[siblingCoElution] = best;
                row.Features[siblingDelta] = nearest;
            }
        }
    }

    /// <summary>
    /// Expected relative intensities of the first <paramref name="count"/> isotopes of a peptide of this neutral mass: a
    /// Poisson with lambda = mass / 1800, close to averagine over tryptic peptide masses.
    /// </summary>
    internal static double[] ExpectedIsotopes(double neutralMass, int count)
    {
        double lambda = neutralMass / 1800;
        var p = new double[count];
        p[0] = Math.Exp(-lambda);
        for (int k = 1; k < count; k++)
            p[k] = p[k - 1] * lambda / k;
        return p;
    }

    /// <summary>
    /// Cosine of an observed MS1 envelope, <paramref name="observedMinus1ToM3"/> (M-1, M0, ... M3), against
    /// <paramref name="expectedM0ToM3"/> with nothing expected at M-1: a peak one isotope below lowers it, as it marks our M0
    /// as probably another ion's isotope. 0 when nothing was observed.
    /// </summary>
    internal static double PeakEnvelopeCosine(double[] observedMinus1ToM3, double[] expectedM0ToM3)
    {
        if (!observedMinus1ToM3.Any(v => v > 0))
            return 0;
        return SpectralSimilarity.CosineOfAlignedVectors(observedMinus1ToM3, [0, .. expectedM0ToM3]);
    }

    private int NearestMs1(double retentionTime)
    {
        int i = Array.BinarySearch(_ms1Rts, retentionTime);
        if (i >= 0)
            return i;
        i = ~i;
        if (i == 0)
            return 0;
        if (i == _ms1Rts.Length)
            return _ms1Rts.Length - 1;
        return retentionTime - _ms1Rts[i - 1] <= _ms1Rts[i] - retentionTime ? i - 1 : i;
    }

    /// <summary>Per-scan median of the traces, each normalised to sum 1 over [from, to]; zero outside the range.</summary>
    private static double[] MedianProfile(IReadOnlyList<double[]> traces, int from, int to, int length)
    {
        var median = new double[length];
        var sums = traces.Select(t => Enumerable.Range(from, to - from + 1).Sum(s => t[s])).ToArray();
        var values = new double[traces.Count];
        for (int s = from; s <= to; s++)
        {
            for (int f = 0; f < traces.Count; f++)
                values[f] = sums[f] > 0 ? traces[f][s] / sums[f] : 0;
            var sorted = values.Order().ToArray();
            median[s] = sorted.Length % 2 == 1 ? sorted[sorted.Length / 2] : 0.5 * (sorted[sorted.Length / 2 - 1] + sorted[sorted.Length / 2]);
        }
        return median;
    }

    /// <summary>Trapezoid area of the trace over [start, end], less the baseline at its lower end point.</summary>
    private static double Area(double[] trace, int start, int end)
    {
        double area = 0;
        for (int s = start; s < end; s++)
            area += 0.5 * (trace[s] + trace[s + 1]);
        double baseline = Math.Min(trace[start], trace[end]) * (end - start);
        return Math.Max(0, area - baseline);
    }

    /// <summary>
    /// Combines each match's features into one score with mzLib's semi-supervised <see cref="TargetDecoyRescorer"/>, with folds
    /// grouped by sequence so that no match is scored by a model trained on it. Each precursor then keeps its best-scoring
    /// candidate, and (D+1)/T q-values are assigned among targets.
    /// If rescoring cannot run (for example too few matches), the pre-rescoring score stands.
    /// </summary>
    private List<DiaPrecursorMatch> AssignQValues(List<DiaPrecursorMatch> matches, (double Low, double High)[] windowBounds, MsDataScan[][] windowScans, double cycleMinutes,
        out List<DiaPrecursorMatch> removed, out Dictionary<int, List<(MzLibUtil.RtMinutes ApexRt, double Score)>> losing)
    {
        var clock = System.Diagnostics.Stopwatch.StartNew();
        // The scores stay in an array beside the candidates until one per precursor is chosen: a copy of every candidate
        // to carry its score would hold two of each, tens of millions, at once
        double[] scores = matches.Select(m => m.Score).ToArray();
        if (matches.Count > 0 && Rescore(matches) is { } rescored)
            scores = rescored;
        Status($"timing: classifier ({_parameters.ClassifierModel}) trained and scored {matches.Count} rows  [{clock.Elapsed:mm\\:ss\\.f}]");

        // One match per precursor: its best-scoring candidate apex, across windows; the others are kept as apex and score
        var candidates = matches;
        var byPrecursor = Enumerable.Range(0, candidates.Count).GroupBy(i => candidates[i].PrecursorIndex)
            .Select(g => g.OrderByDescending(i => scores[i]).ThenBy(i => candidates[i].ApexRt.Value).ToList()).ToList();
        losing = byPrecursor.Where(g => g.Count > 1)
            .ToDictionary(g => candidates[g[0]].PrecursorIndex, g => g.Skip(1).Select(i => (candidates[i].ApexRt, scores[i])).ToList());
        matches = byPrecursor.Select(g => candidates[g[0]] with { Score = scores[g[0]] }).ToList();
        // The losing candidates and their features go now (losing keeps each one's apex and score), so the stages below
        // hold only the winners. The list is the caller's, which reads it no further.
        candidates.Clear();
        candidates.TrimExcess();

        // Interference removal, after DIA-NN: a match a better co-eluting match explains is dropped, target or decoy alike
        removed = [];
        if (_parameters.InterferenceExplainedFragments > 0)
        {
            int WindowOf(double mz)
            {
                for (int w = 0; w < windowBounds.Length; w++)
                    if (mz >= windowBounds[w].Low && mz < windowBounds[w].High)
                        return w;
                return -10;
            }
            float[] FragmentMzs(int index) => _library.GetEntry(index)?.MatchedFragmentIons.OrderByDescending(LibraryIntensity)
                .Take(_parameters.TopFragmentCount).Select(f => f.Mz).ToArray() ?? [];
            double OwnIons(DiaPrecursorMatch winner, DiaPrecursorMatch rival)
            {
                var winnerEntry = _library.GetEntry(winner.PrecursorIndex);
                var rivalEntry = _library.GetEntry(rival.PrecursorIndex);
                if (winnerEntry is null || rivalEntry is null)
                    return double.NaN;
                var winnerIons = ScoredIons(winnerEntry);
                var rivalAllMzs = rivalEntry.MatchedFragmentIons.Select(f => (double)f.Mz).ToList();
                // A winner with much of its own is not judged, however those ions read
                if (DiaInterferenceRemoval.OwnIonShare(winnerIons, rivalAllMzs, _parameters.FragmentTolerancePpm) > _parameters.InterferenceOwnIonShareCeiling)
                    return double.PositiveInfinity;
                var holding = Enumerable.Range(0, windowBounds.Length).Where(w => rival.PrecursorMz >= windowBounds[w].Low && rival.PrecursorMz < windowBounds[w].High).Select(w => windowScans[w]).ToList();
                return DiaInterferenceRemoval.OwnIonCoElution(winnerIons, winnerEntry.MatchedFragmentIons.Select(f => (double)f.Mz).ToList(),
                    ScoredIons(rivalEntry), rivalAllMzs, holding, (winner.ApexRt.Value + rival.ApexRt.Value) / 2,
                    _parameters.FragmentTolerancePpm, _parameters.ApexHalfWidthScans, _parameters.MostIntenseFragmentPeak);
            }
            var kept = DiaInterferenceRemoval.Remove(matches, FragmentMzs, WindowOf, (_parameters.InterferenceRtToleranceCycles ?? _parameters.ApexHalfWidthScans + 1) * cycleMinutes,
                _parameters.FragmentTolerancePpm, _parameters.InterferenceExplainedFragments, Status, _parameters.InterferenceSameMzOnly, OwnIons, _parameters.InterferenceOwnIonFloor);
            var keptIndices = kept.Select(m => m.PrecursorIndex).ToHashSet();
            removed = matches.Where(m => !keptIndices.Contains(m.PrecursorIndex)).ToList();
            matches = kept;
        }

        if (_parameters.ClassifierTreeBlend && matches.Count > 0)
        {
            var blendClock = System.Diagnostics.Stopwatch.StartNew();
            double[] blended = DiaTreeBlend.Blend(matches.Select(m => m.Features).ToList(), matches.Select(m => m.Score).ToList(),
                matches.Select(m => m.IsDecoy).ToList(), matches.Select(m => m.FullSequence).ToList(),
                _parameters.MaxNetworkTrainingRows ?? matches.Count, _parameters.ClassifierSeed, CommonParameters.MaxThreadsToUsePerFile);
            matches = matches.Select((m, i) => m with { Score = blended[i] }).ToList();
            Status($"timing: tree blend over {matches.Count} precursors  [{blendClock.Elapsed:mm\\:ss\\.f}]");
        }

        if (_parameters.NeighbourEvidence && matches.Count > 0)
        {
            var neighbourClock = System.Diagnostics.Stopwatch.StartNew();
            DiaNeighbourEvidence.PrecursorIons? IonsOf(int index) => _library.GetEntry(index) is { } entry
                ? new DiaNeighbourEvidence.PrecursorIons(ScoredIons(entry), entry.MatchedFragmentIons.Select(f => (double)f.Mz).ToList(), entry.BaseSequence)
                : null;
            matches = DiaNeighbourEvidence.Apply(matches, IonsOf, windowBounds, windowScans, (_parameters.InterferenceRtToleranceCycles ?? _parameters.ApexHalfWidthScans + 1) * cycleMinutes,
                _parameters.FragmentTolerancePpm, _parameters.ApexHalfWidthScans, _parameters.MostIntenseFragmentPeak, _parameters.TopFragmentCount, _parameters.ClassifierSeed, CommonParameters.MaxThreadsToUsePerFile, Status);
            Status($"timing: neighbour evidence over {matches.Count} precursors  [{neighbourClock.Elapsed:mm\\:ss\\.f}]");
        }

        var targets = matches.Where(m => !m.IsDecoy).ToList();
        var decoyScores = matches.Where(m => m.IsDecoy).Select(m => m.Score).ToList();
        double[] qValues = DeconvolutionQValueCalculator.AssignQValues(targets.Select(m => m.Score).ToList(), decoyScores);

        // mzLib's calculator gives every target 0 when it has no decoy scores. (D+1)/T still holds here: D is 0, and 1/rank,
        // made monotone, is 1/T for every target.
        if (decoyScores.Count == 0 && targets.Count > 0)
        {
            Warn($"No decoy precursor was scored, so all {targets.Count} target q-values are (0+1)/{targets.Count}, with no decoy evidence behind them. " +
                "Check that the library's decoys have fragments and cover this run's isolation windows and iRT range.");
            Array.Fill(qValues, 1.0 / targets.Count);
        }

        var assigned = targets.Select((m, i) => m with { QValue = qValues[i] }).ToList();
        assigned.AddRange(matches.Where(m => m.IsDecoy));
        return assigned;
    }

    /// <summary>
    /// The classifier's score for every candidate, or null when it could not give them all a finite one. A method of its own,
    /// never inlined, so that the lists it hands the rescorer (one entry per candidate, its features among them) end with its
    /// frame: left in the caller's, they kept every losing candidate's features alive through the stages after.
    /// </summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.NoInlining)]
    private double[]? Rescore(List<DiaPrecursorMatch> matches)
    {
        var rescored = TargetDecoyRescorer.Score(
            matches.Select(m => m.Features).ToList(),
            matches.Select(m => m.IsDecoy).ToList(),
            matches.Select(m => m.FullSequence).ToList(),
            iterations: _parameters.ClassifierLinearIterations,
            positiveQValue: _parameters.ClassifierTrainingQValue,
            model: _parameters.ClassifierModel,
            maxNetworkTrainingRows: _parameters.MaxNetworkTrainingRows,
            randomSeed: _parameters.ClassifierSeed,
            networkMembers: _parameters.ClassifierNetworkMembers,
            networkEpochs: _parameters.ClassifierNetworkEpochs,
            networkPasses: _parameters.ClassifierNetworkPasses,
            networkTrainingSample: _parameters.ClassifierNetworkTrainingSample,
            normalizationGroups: _parameters.ClassifierNormalizationGroups,
            networkLayers: _parameters.EffectiveClassifierNetworkLayers,
            networkPositiveQValue: _parameters.ClassifierNetworkPositiveQValue,
            networkTargetFraction: _parameters.ClassifierNetworkTargetFraction,
            candidateGroups: matches.Select(m => m.PrecursorIndex).ToList());
        return rescored.Scores.All(double.IsFinite) ? rescored.Scores : null;
    }

    /// <summary>The fragments a precursor is scored on, most intense first: the fragment rule's if at least three pass it, as the search picks them.</summary>
    private List<(double Mz, double Intensity)> ScoredIons(MslLibraryEntry entry)
    {
        var ions = entry.MatchedFragmentIons;
        var scorable = _parameters.DiaNnFragmentFilter
            ? ions.Where(f => IsDiaNnScorable(f, _parameters.FragmentMinimumMz, _parameters.FragmentMinimumResidues, _parameters.FragmentMaximumCharge, entry.ChargeState, _parameters.FragmentChargeBelowPrecursor)).ToList()
            : ions.ToList();
        return (scorable.Count >= 3 ? scorable : ions.ToList()).OrderByDescending(LibraryIntensity)
            .Take(_parameters.TopFragmentCount + _parameters.ExtraFragmentCount).Select(f => ((double)f.Mz, LibraryIntensity(f))).ToList();
    }

    /// <summary>A fragment's library intensity, as this run reads it (<see cref="DiaLibrarySearchParameters.FragmentIntensities"/>).</summary>
    private double LibraryIntensity(MslFragmentIon fragment) => _parameters.FragmentIntensities?.Intensity(fragment) ?? fragment.Intensity;
}

/// <summary>
/// The precursors a <see cref="DiaLibrarySearchEngine"/> scored, targets with their q-values, and the peptides they
/// collapse to, with peptide-level q-values from <see cref="DiaPeptideFdr"/>.
/// </summary>
public class DiaLibrarySearchResults(DiaLibrarySearchEngine engine, List<DiaPrecursorMatch> matches, List<DiaPeptideMatch> peptides)
    : MetaMorpheusEngineResults(engine)
{
    public List<DiaPrecursorMatch> Matches { get; init; } = matches;

    public List<DiaPeptideMatch> Peptides { get; init; } = peptides;

    /// <summary>
    /// Matches dropped by interference removal (<see cref="DiaInterferenceRemoval"/>): explained by a better co-eluting
    /// match, so not reported and without a q-value. Kept so that a miss can be told from a precursor never scored.
    /// </summary>
    public List<DiaPrecursorMatch> RemovedAsInterference { get; init; } = [];

    /// <summary>
    /// Each precursor's other candidate peaks (<see cref="DiaLibrarySearchParameters.MaxApexCandidates"/>), as apex RT and
    /// classifier score, best first; precursors with a single candidate are absent. Kept so that a miss where another
    /// engine chose a different peak can be told apart: scored lower among our candidates, or never a candidate. The tree
    /// blend reorders only the chosen candidates, afterwards, so these scores are not comparable with a blended match's.
    /// </summary>
    public IReadOnlyDictionary<int, List<(MzLibUtil.RtMinutes ApexRt, double Score)>> LosingCandidates { get; init; } =
        new Dictionary<int, List<(MzLibUtil.RtMinutes ApexRt, double Score)>>();

    public int TargetCount => Matches.Count(m => !m.IsDecoy);

    public int DecoyCount => Matches.Count(m => m.IsDecoy);

    public override string ToString()
    {
        var sb = new StringBuilder();
        sb.AppendLine(base.ToString());
        sb.AppendLine($"Target precursors: {TargetCount}");
        sb.AppendLine($"Decoy precursors: {DecoyCount}");
        sb.AppendLine($"Target precursors with q-value <= 0.01: {Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01)}");
        sb.AppendLine($"Target peptides with q-value <= 0.01: {Peptides.Count(p => !p.IsDecoy && p.QValue <= 0.01)}");
        return sb.ToString();
    }
}
