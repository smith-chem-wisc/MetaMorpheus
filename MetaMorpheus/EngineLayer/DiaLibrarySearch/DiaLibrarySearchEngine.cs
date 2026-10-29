#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using MzLibUtil;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;

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
    private readonly IIrtMap _irtMap;
    private readonly DiaLibrarySearchParameters _parameters;

    /// <param name="scans">The run's scans. Only MS2 scans with an isolation range are searched.</param>
    /// <param name="library">A loaded library holding targets and decoys. The engine does not dispose it.</param>
    /// <param name="irtMap">This run's retention time to library iRT map.</param>
    public DiaLibrarySearchEngine(MsDataScan[] scans, MslLibrary library, IIrtMap irtMap, DiaLibrarySearchParameters parameters,
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
        var windows = _scans
            .Where(scan => scan.MsnOrder == 2 && scan.IsolationRange is not null)
            .GroupBy(scan => (scan.IsolationRange.Minimum, scan.IsolationRange.Maximum))
            .OrderBy(window => window.Key.Minimum)
            .ToList();

        var best = new Dictionary<int, DiaPrecursorMatch>();
        for (int w = 0; w < windows.Count; w++)
        {
            if (GlobalVariables.StopLoops)
                break;
            var scans = windows[w].OrderBy(scan => scan.RetentionTime).ToArray();
            foreach (var match in SearchWindow(windows[w].Key, scans, tolerance))
                if (!best.TryGetValue(match.PrecursorIndex, out var existing) || match.Score > existing.Score)
                    best[match.PrecursorIndex] = match;
            ReportProgress(new ProgressEventArgs((int)(100.0 * (w + 1) / windows.Count), "Searching DIA windows...", NestedIds));
        }

        var matches = AssignQValues(best.Values.ToList());
        Status("Done.");
        return new DiaLibrarySearchResults(this, matches);
    }

    private List<DiaPrecursorMatch> SearchWindow((double Minimum, double Maximum) window, MsDataScan[] scans, PpmTolerance tolerance)
    {
        var matches = new List<DiaPrecursorMatch>();
        if (scans.Length == 0)
            return matches;

        double[] scanIrts = scans.Select(scan => _irtMap.ToIrt(new RtMinutes(scan.RetentionTime)).Value).ToArray();
        double irtLow = scanIrts.Min() - _parameters.IrtHalfWindow;
        double irtHigh = scanIrts.Max() + _parameters.IrtHalfWindow;

        MslPrecursorIndexEntry[] candidates;
        using (var hits = _library.QueryWindow((float)window.Minimum, (float)window.Maximum, (float)irtLow, (float)irtHigh, includeDecoys: true))
            candidates = hits.Entries.ToArray();

        foreach (var candidate in candidates)
        {
            var entry = _library.GetEntry(candidate.PrecursorIdx);
            if (entry is null || entry.MatchedFragmentIons.Count == 0)
                continue;
            var match = ScoreCandidate(candidate, entry, scans, scanIrts, tolerance);
            if (match is not null)
                matches.Add(match);
        }
        return matches;
    }

    private DiaPrecursorMatch? ScoreCandidate(MslPrecursorIndexEntry candidate, MslLibraryEntry entry, MsDataScan[] scans,
        double[] scanIrts, PpmTolerance tolerance)
    {
        double[] libraryIntensities = entry.MatchedFragmentIons.Select(f => (double)f.Intensity).ToArray();
        double[]? apexIntensities = null;
        double apexSum = 0;
        int apex = -1;

        for (int s = 0; s < scans.Length; s++)
        {
            if (Math.Abs(scanIrts[s] - candidate.Irt) > _parameters.IrtHalfWindow)
                continue;
            double[] observed = FragmentIntensities(scans[s].MassSpectrum, entry.MatchedFragmentIons, tolerance);
            double sum = observed.Sum();
            if (sum > apexSum)
            {
                apexSum = sum;
                apexIntensities = observed;
                apex = s;
            }
        }

        if (apexIntensities is null)
            return null;

        var apexRt = new RtMinutes(scans[apex].RetentionTime);
        return new DiaPrecursorMatch(
            candidate.PrecursorIdx,
            entry.FullSequence,
            candidate.Charge,
            candidate.PrecursorMz,
            candidate.IsDecoy != 0,
            new Irt(candidate.Irt),
            apexRt,
            _irtMap.ToIrt(apexRt),
            SpectralSimilarity.CosineOfAlignedVectors(apexIntensities, libraryIntensities));
    }

    /// <summary>The intensity of the closest peak to each fragment, or 0 when none is within tolerance.</summary>
    private static double[] FragmentIntensities(MzSpectrum spectrum, List<MslFragmentIon> fragments, PpmTolerance tolerance)
    {
        var intensities = new double[fragments.Count];
        if (spectrum.Size == 0)
            return intensities;
        for (int f = 0; f < fragments.Count; f++)
        {
            int i = spectrum.GetClosestPeakIndex(fragments[f].Mz);
            if (tolerance.Within(spectrum.XArray[i], fragments[f].Mz))
                intensities[f] = spectrum.YArray[i];
        }
        return intensities;
    }

    private static List<DiaPrecursorMatch> AssignQValues(List<DiaPrecursorMatch> matches)
    {
        var targets = matches.Where(m => !m.IsDecoy).ToList();
        var decoyScores = matches.Where(m => m.IsDecoy).Select(m => m.Score).ToList();
        double[] qValues = DeconvolutionQValueCalculator.AssignQValues(targets.Select(m => m.Score).ToList(), decoyScores);

        var assigned = targets.Select((m, i) => m with { QValue = qValues[i] }).ToList();
        assigned.AddRange(matches.Where(m => m.IsDecoy));
        return assigned;
    }
}

/// <summary>The precursors a <see cref="DiaLibrarySearchEngine"/> scored, targets with their q-values.</summary>
public class DiaLibrarySearchResults(DiaLibrarySearchEngine engine, List<DiaPrecursorMatch> matches) : MetaMorpheusEngineResults(engine)
{
    public List<DiaPrecursorMatch> Matches { get; init; } = matches;

    public int TargetCount => Matches.Count(m => !m.IsDecoy);

    public int DecoyCount => Matches.Count(m => m.IsDecoy);

    public override string ToString()
    {
        var sb = new StringBuilder();
        sb.AppendLine(base.ToString());
        sb.AppendLine($"Target precursors: {TargetCount}");
        sb.AppendLine($"Decoy precursors: {DecoyCount}");
        sb.AppendLine($"Target precursors with q-value <= 0.01: {Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01)}");
        return sb.ToString();
    }
}
