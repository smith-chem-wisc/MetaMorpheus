using System;
using System.Collections.Generic;
using System.Linq;
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
    }

    /// <summary>
    /// Pass 2 indexed dual single-series scoring engine (01_Architecture.md decisions #4-#8). Builds a
    /// fragment index from the parent proteoforms (following the ModernSearchEngine binning pattern,
    /// <see cref="FragmentBinsPerDalton"/>), then for each MS2 scan that Pass 1 did not match intact
    /// (#4a) gathers candidate parents via the index and the <see cref="TruncationAcceptor"/>, and scores
    /// each candidate twice — once over its N-terminal-ion series, once over its C-terminal-ion series —
    /// using the standard MetaMorpheus matched-ion score. The highest-scoring (parent, series) pair wins.
    ///
    /// Opposing-series ions are simply never matched in a given direction, so they neither disqualify a
    /// candidate nor contribute to its score (#8). No chopping or FDR here — that is Pass 3.
    /// </summary>
    public class TruncationSearchEngine
    {
        public const int FragmentBinsPerDalton = 1000;

        /// <summary>Existing SearchParameters.DefaultMaxFragmentSize (decision #6).</summary>
        public const double DefaultMaxFragmentSize = 30000;

        private readonly List<TruncationParent> _parents; // included parents, sorted ascending by mass
        private readonly Ms2ScanWithSpecificMass[] _scans;
        private readonly CommonParameters _commonParameters;
        private readonly TruncationAcceptor _acceptor;
        private readonly IReadOnlyDictionary<int, double> _pass1TheoreticalMassByScanNumber;
        private List<int>[] _fragmentIndex;

        /// <summary>Non-fatal anomalies surfaced to the caller (e.g. the MaxFragmentSize exclusion summary, #6).</summary>
        public List<string> Warnings { get; } = new();

        /// <summary>Number of parents excluded because their mass exceeded MaxFragmentSize (decision #6).</summary>
        public int ExcludedOversizedParentCount { get; }

        public TruncationSearchEngine(
            IEnumerable<TruncationParent> parents,
            Ms2ScanWithSpecificMass[] scans,
            CommonParameters commonParameters,
            TruncationAcceptor acceptor,
            IReadOnlyDictionary<int, double> pass1TheoreticalMassByScanNumber = null,
            double maxFragmentSize = DefaultMaxFragmentSize)
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
        }

        public List<TruncationParentSelection> Run()
        {
            BuildFragmentIndex();

            var results = new List<TruncationParentSelection>(_scans.Length);
            var nTermProducts = new List<Product>();
            var cTermProducts = new List<Product>();

            for (int scanIndex = 0; scanIndex < _scans.Length; scanIndex++)
            {
                Ms2ScanWithSpecificMass scan = _scans[scanIndex];

                // Intact-skip (#4a): Pass 1 matched this scan and the parent's theoretical mass equals
                // the scan precursor mass within tolerance. No Pass 2 scoring; carry it through as intact.
                if (_pass1TheoreticalMassByScanNumber != null
                    && _pass1TheoreticalMassByScanNumber.TryGetValue(scan.OneBasedScanNumber, out double pass1Mass)
                    && _commonParameters.PrecursorMassTolerance.Within(scan.PrecursorMass, pass1Mass))
                {
                    results.Add(new TruncationParentSelection
                    {
                        ScanIndex = scanIndex,
                        Scan = scan,
                        Outcome = TruncationScanOutcome.IntactInherited
                    });
                    continue;
                }

                List<int> candidateParentIds = GatherCandidateParents(scan);
                if (candidateParentIds.Count == 0)
                {
                    results.Add(NoWinner(scanIndex, scan));
                    continue;
                }

                TruncationParent bestParent = null;
                FragmentationTerminus bestSeries = FragmentationTerminus.None;
                double bestScore = 0;

                foreach (int id in candidateParentIds)
                {
                    TruncationParent parent = _parents[id];

                    parent.Proteoform.Fragment(_commonParameters.DissociationType, FragmentationTerminus.N, nTermProducts, _commonParameters.FragmentationParameters);
                    double nScore = MetaMorpheusEngine.CalculatePeptideScore(scan.TheScan,
                        MetaMorpheusEngine.MatchFragmentIons(scan, nTermProducts, _commonParameters));

                    parent.Proteoform.Fragment(_commonParameters.DissociationType, FragmentationTerminus.C, cTermProducts, _commonParameters.FragmentationParameters);
                    double cScore = MetaMorpheusEngine.CalculatePeptideScore(scan.TheScan,
                        MetaMorpheusEngine.MatchFragmentIons(scan, cTermProducts, _commonParameters));

                    if (nScore > bestScore)
                    {
                        bestScore = nScore;
                        bestParent = parent;
                        bestSeries = FragmentationTerminus.N;
                    }

                    if (cScore > bestScore)
                    {
                        bestScore = cScore;
                        bestParent = parent;
                        bestSeries = FragmentationTerminus.C;
                    }
                }

                if (bestParent != null && bestScore >= _commonParameters.ScoreCutoff)
                {
                    results.Add(new TruncationParentSelection
                    {
                        ScanIndex = scanIndex,
                        Scan = scan,
                        Outcome = TruncationScanOutcome.Winner,
                        WinningParent = bestParent,
                        WinningSeries = bestSeries,
                        Score = bestScore
                    });
                }
                else
                {
                    results.Add(NoWinner(scanIndex, scan));
                }
            }

            return results;
        }

        private static TruncationParentSelection NoWinner(int scanIndex, Ms2ScanWithSpecificMass scan) =>
            new() { ScanIndex = scanIndex, Scan = scan, Outcome = TruncationScanOutcome.NoWinner };

        /// <summary>
        /// Builds the fragment index over both ion series of every included parent, mapping each
        /// rounded fragment-mass bin to the ids of parents that have a fragment in that bin. Parents are
        /// stored in ascending-mass order so each bin's id list is mass-sorted.
        /// </summary>
        private void BuildFragmentIndex()
        {
            if (_parents.Count == 0)
            {
                _fragmentIndex = Array.Empty<List<int>>();
                return;
            }

            double maxParentMass = _parents[_parents.Count - 1].MonoisotopicMass;
            int indexLength = (int)Math.Ceiling(maxParentMass * FragmentBinsPerDalton) + 2;
            _fragmentIndex = new List<int>[indexLength];

            var products = new List<Product>();
            for (int id = 0; id < _parents.Count; id++)
            {
                _parents[id].Proteoform.Fragment(_commonParameters.DissociationType, FragmentationTerminus.Both, products, _commonParameters.FragmentationParameters);
                foreach (Product product in products)
                {
                    if (double.IsNaN(product.NeutralMass))
                    {
                        continue;
                    }

                    int bin = (int)Math.Round(product.NeutralMass * FragmentBinsPerDalton);
                    if (bin < 0 || bin >= _fragmentIndex.Length)
                    {
                        continue;
                    }

                    (_fragmentIndex[bin] ??= new List<int>()).Add(id);
                }
            }
        }

        /// <summary>
        /// Gathers candidate parent ids for a scan: any parent that has a fragment within product
        /// tolerance of an experimental fragment AND passes the precursor acceptor (#7). The acceptor
        /// excludes parents lighter than the scan precursor — they cannot be truncated up to it.
        /// </summary>
        private List<int> GatherCandidateParents(Ms2ScanWithSpecificMass scan)
        {
            var candidateIds = new HashSet<int>();

            if (_fragmentIndex.Length == 0 || scan.ExperimentalFragments == null)
            {
                return candidateIds.ToList();
            }

            foreach (IsotopicEnvelope envelope in scan.ExperimentalFragments)
            {
                double experimentalMass = envelope.MonoisotopicMass;

                int floorBin = Math.Max(0,
                    (int)Math.Floor(_commonParameters.ProductMassTolerance.GetMinimumValue(experimentalMass) * FragmentBinsPerDalton));
                int ceilingBin = Math.Min(_fragmentIndex.Length - 1,
                    (int)Math.Ceiling(_commonParameters.ProductMassTolerance.GetMaximumValue(experimentalMass) * FragmentBinsPerDalton));

                for (int b = floorBin; b <= ceilingBin; b++)
                {
                    List<int> bin = _fragmentIndex[b];
                    if (bin == null)
                    {
                        continue;
                    }

                    foreach (int id in bin)
                    {
                        if (!candidateIds.Contains(id) && _acceptor.Accepts(scan.PrecursorMass, _parents[id].MonoisotopicMass) >= 0)
                        {
                            candidateIds.Add(id);
                        }
                    }
                }
            }

            return candidateIds.ToList();
        }
    }
}
