using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MassSpectrometry;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer.Truncation
{
    /// <summary>What happened for a given scan in Pass 2.</summary>
    public enum TruncationScanOutcome
    {
        /// <summary>Pass 1 already matched this scan intact (decision #4a); no Pass 2 scoring was run.</summary>
        IntactInherited,

        /// <summary>A (parent, terminus) pair won the dual single-series scoring above the score cutoff.</summary>
        Winner,

        /// <summary>No parent passed the precursor acceptor and score cutoff for this scan.</summary>
        NoWinner
    }

    /// <summary>
    /// A full-length parent proteoform sourced from Pass 1, carrying the provenance needed to tie
    /// pipe-ambiguous alternatives back to one Pass 1 row (decision #2) and to report protein accession.
    /// </summary>
    public class TruncationParent
    {
        public TruncationParent(PeptideWithSetModifications proteoform, string proteinAccession, object originatingId, bool isDecoy)
        {
            Proteoform = proteoform;
            ProteinAccession = proteinAccession;
            OriginatingId = originatingId;
            IsDecoy = isDecoy;
        }

        public PeptideWithSetModifications Proteoform { get; }
        public string ProteinAccession { get; }
        public object OriginatingId { get; }
        public bool IsDecoy { get; }
        public double MonoisotopicMass => Proteoform.MonoisotopicMass;
    }

    /// <summary>
    /// Result of Pass 2 for a single scan: the winning parent and which ion series matched it, or an
    /// intact-inherited / no-winner outcome. The chopped terminus (Pass 3) is the OPPOSITE of
    /// <see cref="WinningSeries"/> (a winning N-terminal-ion series implies the C-terminus was lost).
    /// </summary>
    public class TruncationParentSelection
    {
        public int ScanIndex { get; init; }
        public Ms2ScanWithSpecificMass Scan { get; init; }
        public TruncationScanOutcome Outcome { get; init; }
        public TruncationParent WinningParent { get; init; }
        public FragmentationTerminus WinningSeries { get; init; } = FragmentationTerminus.None;
        public double Score { get; init; }

        /// <summary>0-based rank of the winning parent in the per-scan index match-count ordering (diagnostic:
        /// how deep into the rescored top-N the winner sat). -1 when not a Winner.</summary>
        public int CandidateRank { get; init; } = -1;

        /// <summary>How many candidate parents passed the match-count prune for this scan (diagnostic).</summary>
        public int CandidatePoolSize { get; init; }
    }

    /// <summary>
    /// Pass 2 indexed dual single-series scoring engine (01_Architecture.md decisions #4-#8). Builds two
    /// fragment indexes — one per ion series (N: b/c, C: y/z) — from the parent proteoforms
    /// (<see cref="FragmentBinsPerDalton"/> binning, the ModernSearchEngine pattern), then for each MS2 scan
    /// that Pass 1 did not match intact (#4a) it scores candidate parents directly off the index: a single
    /// pass over the scan's experimental peaks accumulates, per parent, how many of its N-series and
    /// C-series fragment bins are hit. This scales to a whole-database parent set because no parent is
    /// re-fragmented during candidate gathering. Parents whose best single-series match count is too low to
    /// possibly clear <see cref="CommonParameters.ScoreCutoff"/> are pruned; the top remaining candidates are
    /// then authoritatively re-scored with the standard MetaMorpheus matched-ion score (over each series
    /// separately, #8/#12), and the highest-scoring (parent, series) pair wins. No chopping or FDR here — Pass 3.
    /// </summary>
    public class TruncationSearchEngine
    {
        public const int FragmentBinsPerDalton = 1000;

        /// <summary>Existing SearchParameters.DefaultMaxFragmentSize (decision #6).</summary>
        public const double DefaultMaxFragmentSize = 30000;

        /// <summary>
        /// Cap on how many index-ranked candidates per scan get the authoritative full both-series rescoring.
        /// The true highest-scoring parent has the most matched fragments, so it sits at the top of the
        /// match-count ranking; rescoring a generous slice makes the (parent, series) winner identical to the
        /// brute-force result in all but pathological cases while bounding per-scan work.
        /// </summary>
        public const int MaxCandidatesToScore = 50;

        private readonly List<TruncationParent> _parents; // included parents, sorted ascending by mass
        private readonly Ms2ScanWithSpecificMass[] _scans;
        private readonly CommonParameters _commonParameters;
        private readonly TruncationAcceptor _acceptor;
        private readonly IReadOnlyDictionary<int, double> _pass1TheoreticalMassByScanNumber;
        private readonly DissociationType[] _dissociationTypes; // index is built over all of these (#5)
        private readonly int _minMatchedFragments; // below this a single series cannot reach ScoreCutoff
        private List<int>[] _nIndex; // N-series (b/c) fragment-mass bin -> ids of parents with a fragment there
        private List<int>[] _cIndex; // C-series (y/z) fragment-mass bin -> ids of parents with a fragment there

        // Optional per-scan candidate restriction (tag-filter per-scan mode): a scan is scored only against
        // parents whose underlying (DECOY-stripped) protein accession is in the scan's allowed set, so a scan
        // cannot be won by a protein its own sequence tags did not support. Null = no restriction (default).
        private readonly IReadOnlyDictionary<Ms2ScanWithSpecificMass, HashSet<string>> _allowedAccessionsByScan;
        private static readonly HashSet<string> NoAllowedProteins = new();

        /// <summary>Non-fatal anomalies surfaced to the caller (e.g. the MaxFragmentSize exclusion summary, #6).</summary>
        public List<string> Warnings { get; } = new();

        /// <summary>Number of parents excluded because their mass exceeded MaxFragmentSize (decision #6).</summary>
        public int ExcludedOversizedParentCount { get; }

        /// <summary>Wall-clock seconds spent building the fragment index (perf logging, 03_Benchmarks).</summary>
        public double IndexBuildSeconds { get; private set; }

        /// <summary>Wall-clock seconds spent in the per-scan dual single-series scoring loop.</summary>
        public double ScoringSeconds { get; private set; }

        /// <summary>Count of parents actually placed in the index (after the oversize exclusion).</summary>
        public int IndexedParentCount => _parents.Count;

        public TruncationSearchEngine(
            IEnumerable<TruncationParent> parents,
            Ms2ScanWithSpecificMass[] scans,
            CommonParameters commonParameters,
            TruncationAcceptor acceptor,
            IReadOnlyDictionary<int, double> pass1TheoreticalMassByScanNumber = null,
            double maxFragmentSize = DefaultMaxFragmentSize,
            IReadOnlyDictionary<Ms2ScanWithSpecificMass, HashSet<string>> allowedProteinAccessionsByScan = null)
        {
            var allParents = parents.ToList();

            // Exclude oversized parents from the index; emit a single summary line (#6).
            _parents = allParents
                .Where(p => p.MonoisotopicMass <= maxFragmentSize)
                .OrderBy(p => p.MonoisotopicMass)
                .ToList();
            ExcludedOversizedParentCount = allParents.Count - _parents.Count;
            if (ExcludedOversizedParentCount > 0)
            {
                Warnings.Add($"{ExcludedOversizedParentCount} parent proteoforms exceeded MaxFragmentSize and were excluded.");
            }

            _scans = scans;
            _commonParameters = commonParameters;
            _acceptor = acceptor;
            _pass1TheoreticalMassByScanNumber = pass1TheoreticalMassByScanNumber;
            _allowedAccessionsByScan = allowedProteinAccessionsByScan;

            // A single series scores 1 + intensityFraction (<2) per matched ion, so a parent matching fewer
            // than floor(ScoreCutoff/2)+1 fragments in a series can never reach ScoreCutoff and is safe to
            // prune before the authoritative rescoring (result-preserving for the default cutoff).
            _minMatchedFragments = Math.Max(1, (int)Math.Floor(_commonParameters.ScoreCutoff / 2.0) + 1);

            // When the dissociation type is Autodetect it is resolved per scan from the scan header
            // (mirroring ClassicSearchEngine). The index must hold fragments for every type the scans use.
            _dissociationTypes = _commonParameters.DissociationType != DissociationType.Autodetect
                ? new[] { _commonParameters.DissociationType }
                : _scans.Select(ResolveDissociationType).Distinct().ToArray();
        }

        /// <summary>
        /// Resolves the effective dissociation type for a scan: the fixed CommonParameters type, or — when
        /// that is Autodetect — the scan-header type (falling back to HCD when the header omits it).
        /// </summary>
        public static DissociationType ResolveDissociationType(CommonParameters commonParameters, Ms2ScanWithSpecificMass scan) =>
            commonParameters.DissociationType != DissociationType.Autodetect
                ? commonParameters.DissociationType
                : scan.TheScan.DissociationType ?? DissociationType.HCD;

        private DissociationType ResolveDissociationType(Ms2ScanWithSpecificMass scan) =>
            ResolveDissociationType(_commonParameters, scan);

        public List<TruncationParentSelection> Run()
        {
            var stopwatch = System.Diagnostics.Stopwatch.StartNew();
            BuildFragmentIndex();
            IndexBuildSeconds = stopwatch.Elapsed.TotalSeconds;
            stopwatch.Restart();

            // Per-scan scoring is independent, so parallelize across scans (results land in a fixed slot to
            // preserve order). The two indexes and the acceptor are read-only after the build, so they are
            // shared safely; each iteration uses its own product/accumulator buffers.
            var results = new TruncationParentSelection[_scans.Length];
            var parallelOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = Math.Max(1, _commonParameters.MaxThreadsToUsePerFile)
            };

            Parallel.For(0, _scans.Length, parallelOptions, scanIndex =>
            {
                Ms2ScanWithSpecificMass scan = _scans[scanIndex];

                // Intact-skip (#4a): Pass 1 matched this scan and the parent's theoretical mass equals the
                // scan precursor mass within tolerance. No Pass 2 scoring; carry it through as intact.
                if (_pass1TheoreticalMassByScanNumber != null
                    && _pass1TheoreticalMassByScanNumber.TryGetValue(scan.OneBasedScanNumber, out double pass1Mass)
                    && _commonParameters.PrecursorMassTolerance.Within(scan.PrecursorMass, pass1Mass))
                {
                    results[scanIndex] = new TruncationParentSelection
                    {
                        ScanIndex = scanIndex,
                        Scan = scan,
                        Outcome = TruncationScanOutcome.IntactInherited
                    };
                    return;
                }

                // Index-based candidate selection (#7): rank parents by indexed single-series match count,
                // pruned to those that could possibly clear the cutoff.
                List<int> candidateIds = RankCandidateParents(scan);
                if (candidateIds.Count == 0)
                {
                    results[scanIndex] = NoWinner(scanIndex, scan);
                    return;
                }

                // Authoritative rescoring of the shortlist with the standard both-series matched-ion score,
                // scored separately per series (Pass 2 selects parent + series, #8/#12).
                var nTermProducts = new List<Product>();
                var cTermProducts = new List<Product>();
                TruncationParent bestParent = null;
                FragmentationTerminus bestSeries = FragmentationTerminus.None;
                double bestScore = 0;
                int bestRank = -1;

                DissociationType dissociationType = ResolveDissociationType(scan);

                for (int rank = 0; rank < candidateIds.Count; rank++)
                {
                    TruncationParent parent = _parents[candidateIds[rank]];

                    parent.Proteoform.Fragment(dissociationType, FragmentationTerminus.N, nTermProducts, _commonParameters.FragmentationParameters);
                    double nScore = MetaMorpheusEngine.CalculatePeptideScore(scan.TheScan,
                        MetaMorpheusEngine.MatchFragmentIons(scan, nTermProducts, _commonParameters));

                    parent.Proteoform.Fragment(dissociationType, FragmentationTerminus.C, cTermProducts, _commonParameters.FragmentationParameters);
                    double cScore = MetaMorpheusEngine.CalculatePeptideScore(scan.TheScan,
                        MetaMorpheusEngine.MatchFragmentIons(scan, cTermProducts, _commonParameters));

                    if (nScore > bestScore)
                    {
                        bestScore = nScore;
                        bestParent = parent;
                        bestSeries = FragmentationTerminus.N;
                        bestRank = rank;
                    }

                    if (cScore > bestScore)
                    {
                        bestScore = cScore;
                        bestParent = parent;
                        bestSeries = FragmentationTerminus.C;
                        bestRank = rank;
                    }
                }

                results[scanIndex] = bestParent != null && bestScore >= _commonParameters.ScoreCutoff
                    ? new TruncationParentSelection
                    {
                        ScanIndex = scanIndex,
                        Scan = scan,
                        Outcome = TruncationScanOutcome.Winner,
                        WinningParent = bestParent,
                        WinningSeries = bestSeries,
                        Score = bestScore,
                        CandidateRank = bestRank,
                        CandidatePoolSize = candidateIds.Count
                    }
                    : NoWinner(scanIndex, scan);
            });

            ScoringSeconds = stopwatch.Elapsed.TotalSeconds;
            return results.ToList();
        }

        private static TruncationParentSelection NoWinner(int scanIndex, Ms2ScanWithSpecificMass scan) =>
            new() { ScanIndex = scanIndex, Scan = scan, Outcome = TruncationScanOutcome.NoWinner };

        /// <summary>
        /// Scans the experimental peaks once and, off the prebuilt indexes, counts how many distinct N-series
        /// and C-series fragment bins of each (precursor-acceptable) parent are hit. Parents whose best
        /// single-series count clears <see cref="_minMatchedFragments"/> are returned, highest count first and
        /// truncated to <see cref="MaxCandidatesToScore"/>, for authoritative rescoring. No re-fragmentation.
        /// </summary>
        private List<int> RankCandidateParents(Ms2ScanWithSpecificMass scan)
        {
            if (_nIndex.Length == 0 || scan.ExperimentalFragments == null)
            {
                return new List<int>();
            }

            var nCount = new Dictionary<int, int>();
            var cCount = new Dictionary<int, int>();
            var countedBins = new HashSet<int>();

            // Per-scan candidate restriction (null = none). An empty set means the scan's tags supported no
            // protein, so it gets no candidates.
            HashSet<string> allowed = _allowedAccessionsByScan == null
                ? null
                : (_allowedAccessionsByScan.TryGetValue(scan, out HashSet<string> s) ? s : NoAllowedProteins);

            foreach (IsotopicEnvelope envelope in scan.ExperimentalFragments)
            {
                double experimentalMass = envelope.MonoisotopicMass;
                int floorBin = Math.Max(0,
                    (int)Math.Floor(_commonParameters.ProductMassTolerance.GetMinimumValue(experimentalMass) * FragmentBinsPerDalton));
                int ceilingBin = (int)Math.Ceiling(_commonParameters.ProductMassTolerance.GetMaximumValue(experimentalMass) * FragmentBinsPerDalton);

                for (int b = floorBin; b <= ceilingBin; b++)
                {
                    if (!countedBins.Add(b))
                    {
                        continue; // each fragment-mass bin contributes at most once per parent per scan
                    }

                    Tally(_nIndex, b, scan.PrecursorMass, nCount, allowed);
                    Tally(_cIndex, b, scan.PrecursorMass, cCount, allowed);
                }
            }

            var ids = new HashSet<int>(nCount.Keys);
            ids.UnionWith(cCount.Keys);

            var ranked = new List<(int Id, int Count)>();
            foreach (int id in ids)
            {
                nCount.TryGetValue(id, out int n);
                cCount.TryGetValue(id, out int c);
                int best = Math.Max(n, c);
                if (best >= _minMatchedFragments)
                {
                    ranked.Add((id, best));
                }
            }

            ranked.Sort((a, b) => b.Count.CompareTo(a.Count));
            int take = Math.Min(ranked.Count, MaxCandidatesToScore);
            var shortlist = new List<int>(take);
            for (int i = 0; i < take; i++)
            {
                shortlist.Add(ranked[i].Id);
            }

            return shortlist;
        }

        /// <summary>Increments the match count for every precursor-acceptable, scan-allowed parent in the bin.</summary>
        private void Tally(List<int>[] index, int bin, double precursorMass, Dictionary<int, int> count, HashSet<string> allowed)
        {
            if (bin < 0 || bin >= index.Length)
            {
                return;
            }

            List<int> binList = index[bin];
            if (binList == null)
            {
                return;
            }

            foreach (int id in binList)
            {
                // The acceptor excludes parents lighter than the scan precursor — they cannot be truncated up to it (#7);
                // the per-scan filter (if any) excludes parents whose protein the scan's tags did not support.
                if (_acceptor.Accepts(precursorMass, _parents[id].MonoisotopicMass) >= 0 && IsScanAllowed(id, allowed))
                {
                    count[id] = (count.TryGetValue(id, out int v) ? v : 0) + 1;
                }
            }
        }

        /// <summary>
        /// Per-scan restriction test: a parent is allowed when there is no restriction (<paramref name="allowed"/>
        /// is null), or its underlying target protein accession (stripping any DECOY_ prefix so a target and its
        /// reverse decoy are allowed together, preserving FDR balance) is in the scan's allowed set.
        /// </summary>
        private bool IsScanAllowed(int id, HashSet<string> allowed)
        {
            if (allowed == null)
            {
                return true;
            }
            string accession = _parents[id].ProteinAccession ?? string.Empty;
            if (accession.StartsWith("DECOY_", StringComparison.Ordinal))
            {
                accession = accession.Substring("DECOY_".Length);
            }
            return allowed.Contains(accession);
        }

        /// <summary>
        /// Builds the two series-tagged fragment indexes (N: b/c, C: y/z) over every included parent: each
        /// rounded fragment-mass bin maps to the ids of parents that have a fragment of that series in it.
        /// Fragmenting parents is independent, so this runs in parallel — parents (sorted ascending by mass,
        /// hence id) are split into contiguous chunks, each chunk builds thread-local sparse bin maps, and the
        /// chunks are merged back in id order so every bin's id list stays ascending.
        /// </summary>
        private void BuildFragmentIndex()
        {
            if (_parents.Count == 0)
            {
                _nIndex = Array.Empty<List<int>>();
                _cIndex = Array.Empty<List<int>>();
                return;
            }

            double maxParentMass = _parents[_parents.Count - 1].MonoisotopicMass;
            int indexLength = (int)Math.Ceiling(maxParentMass * FragmentBinsPerDalton) + 2;
            _nIndex = new List<int>[indexLength];
            _cIndex = new List<int>[indexLength];

            int chunks = Math.Max(1, Math.Min(_commonParameters.MaxThreadsToUsePerFile, _parents.Count));
            var partialsN = new Dictionary<int, List<int>>[chunks];
            var partialsC = new Dictionary<int, List<int>>[chunks];

            Parallel.For(0, chunks, chunk =>
            {
                int lo = (int)((long)_parents.Count * chunk / chunks);
                int hi = (int)((long)_parents.Count * (chunk + 1) / chunks);
                var localN = new Dictionary<int, List<int>>();
                var localC = new Dictionary<int, List<int>>();
                var products = new List<Product>();
                var binsN = new HashSet<int>();
                var binsC = new HashSet<int>();

                for (int id = lo; id < hi; id++)
                {
                    binsN.Clear();
                    binsC.Clear();
                    foreach (DissociationType dissociationType in _dissociationTypes)
                    {
                        _parents[id].Proteoform.Fragment(dissociationType, FragmentationTerminus.Both, products, _commonParameters.FragmentationParameters);
                        foreach (Product product in products)
                        {
                            if (double.IsNaN(product.NeutralMass))
                            {
                                continue;
                            }

                            int bin = (int)Math.Round(product.NeutralMass * FragmentBinsPerDalton);
                            if (bin < 0 || bin >= indexLength)
                            {
                                continue;
                            }

                            if (product.Terminus == FragmentationTerminus.N)
                            {
                                binsN.Add(bin);
                            }
                            else if (product.Terminus == FragmentationTerminus.C)
                            {
                                binsC.Add(bin);
                            }
                        }
                    }

                    AddToLocal(localN, binsN, id);
                    AddToLocal(localC, binsC, id);
                }

                partialsN[chunk] = localN;
                partialsC[chunk] = localC;
            });

            MergePartials(partialsN, _nIndex);
            MergePartials(partialsC, _cIndex);
        }

        private static void AddToLocal(Dictionary<int, List<int>> local, HashSet<int> bins, int id)
        {
            foreach (int bin in bins)
            {
                if (!local.TryGetValue(bin, out List<int> binList))
                {
                    local[bin] = binList = new List<int>();
                }
                binList.Add(id); // ids ascending within a chunk's contiguous range
            }
        }

        /// <summary>Merges thread-local bin maps into the global index in chunk (id) order, keeping bins ascending.</summary>
        private static void MergePartials(Dictionary<int, List<int>>[] partials, List<int>[] index)
        {
            for (int chunk = 0; chunk < partials.Length; chunk++)
            {
                foreach (KeyValuePair<int, List<int>> kv in partials[chunk])
                {
                    List<int> global = index[kv.Key];
                    if (global == null)
                    {
                        index[kv.Key] = kv.Value;
                    }
                    else
                    {
                        global.AddRange(kv.Value);
                    }
                }

                partials[chunk] = null; // release as we go
            }
        }
    }
}
