using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.Util;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;

namespace EngineLayer.ParallelSearch
{

    /// <summary>
    /// Designed to be a more memory efficient version of ClassicSearchEngine for use in parallel searches.
    /// The base PSM list is shared across threads, and only the PSMs that are updated by the transient search are cloned and modified, while the rest of the PSMs remain shared and read-only.
    /// Features that are removed from classic search engine for runtime efficiency and memory efficiency:
    /// - Some detailed logging
    /// - Mass Differance acceptor is always exact. 
    /// </summary>
    public class TransientClassicSearchEngine : ClassicSearchEngine
    {
        private bool _singleThreadMode;
        private readonly ConcurrentDictionary<int, byte> UpdatedIndexes = new();
        private readonly bool _copyOnWriteEnabled;
        private readonly PerFileTolerances _perFileTolerances;
        private readonly bool _filterPrecursorPerFile;
        private readonly double[] _scanPrecursorHalfWidth;
        public TransientClassicSearchEngine(SpectralMatch[] globalPsms, Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans,
            List<Modification> variableModifications, List<Modification> fixedModifications, 
            List<IBioPolymer> proteinList, MassDiffAcceptor searchMode, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIds,
            bool copyOnWriteEnabled = false, PerFileTolerances perFileTolerances = null, ClassicSearchScanIndex scanIndex = null)
            : base(globalPsms, arrayOfSortedMS2Scans, variableModifications, fixedModifications, null, null, null, proteinList, searchMode, commonParameters, fileSpecificParameters, null, nestedIds, false,
                scanIndex: scanIndex)
        {
            UpdatedIndexes = new ConcurrentDictionary<int, byte>();
            _singleThreadMode = CommonParameters.MaxThreadsToUsePerFile <= 1;
            _copyOnWriteEnabled = copyOnWriteEnabled;
            _perFileTolerances = perFileTolerances;
            _filterPrecursorPerFile = perFileTolerances?.PrecursorVaries == true;
            _scanPrecursorHalfWidth = perFileTolerances?.ScanPrecursorHalfWidth;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double proteinsSearched = 0;
            int oldPercentProgress = 0;
            int peptideCounter = 0;

            Status("Performing classic search...");

            if (Proteins.Any())
            {
                // as in ClassicSearchEngine: the scans whose candidates can be scored before their ions are built, decided once
                bool[] scoreBeforeMatching = ScanIndex?.ScoreBeforeMatching ?? new bool[ArrayOfSortedMS2Scans.Length];
                if (ScanIndex == null)
                {
                    for (int s = 0; s < scoreBeforeMatching.Length; s++)
                    {
                        scoreBeforeMatching[s] = CanScoreBeforeMatching(ArrayOfSortedMS2Scans[s], CommonParameters, matchAllCharges: false);
                    }
                }

                Action<int, int> processProteinRange = (start, end) =>
                {
                    var peptideTheorProducts = new List<Product>();
                    var sortedProductMasses = new SortedProductMasses();
                    var fragmentMatchScratch = new FragmentMatchScratch();
                    var acceptableScans = new List<ScanWithIndexAndNotchInfo>();
                    int peptidesInRange = 0;
                    for (int i = start; i < end; i++)
                    {
                        // Stop loop if canceled
                        if (GlobalVariables.StopLoops) { break; }

                        // digest each protein into peptides and search for each peptide in all spectra within precursor mass tolerance
                        foreach (var specificBioPolymer in Proteins[i].Digest(CommonParameters.DigestionParams, FixedModifications, VariableModifications))
                        {
                            peptidesInRange++;

                            // fragments are generated on the first scan in the precursor window; most peptides have none
                            peptideTheorProducts.Clear();
                            bool fragmented = false;

                            FillAcceptableScans(specificBioPolymer.MonoisotopicMass, acceptableScans);
                            foreach (ScanWithIndexAndNotchInfo scan in acceptableScans)
                            {
                                if (!fragmented)
                                {
                                    specificBioPolymer.Fragment(CommonParameters.DissociationType, CommonParameters.DigestionParams.FragmentationTerminus, peptideTheorProducts, CommonParameters.FragmentationParameters);
                                    // sorted once per peptide, so that each scan matches all of them in one walk of its masses
                                    sortedProductMasses.Build(peptideTheorProducts);
                                    fragmented = true;
                                }

                                Ms2ScanWithSpecificMass theScan = ArrayOfSortedMS2Scans[scan.ScanIndex];
                                CommonParameters scoringParameters = _perFileTolerances == null
                                    ? CommonParameters
                                    : _perFileTolerances.Parameters[_perFileTolerances.ScanFileIndex[scan.ScanIndex]];

                                List<MatchedFragmentIon> matchedIons;
                                double score;
                                if (scoreBeforeMatching[scan.ScanIndex])
                                {
                                    // the ions are built only when the score reaches the cutoff (null below it)
                                    score = ScoreFragmentIons(theScan, peptideTheorProducts, sortedProductMasses, scoringParameters,
                                        CommonParameters.ScoreCutoff, fragmentMatchScratch, out matchedIons);
                                }
                                else
                                {
                                    matchedIons = MatchFragmentIons(theScan, peptideTheorProducts, scoringParameters);
                                    score = CalculatePeptideScore(theScan.TheScan, matchedIons);
                                }

                                if (score < CommonParameters.ScoreCutoff)
                                    continue;

                                AddPeptideCandidateToPsm(scan, score, specificBioPolymer, matchedIons);
                            }
                        }

                        // report search progress (proteins searched so far out of total proteins in database)
                        proteinsSearched++;
                        var percentProgress = (int)((proteinsSearched / Proteins.Count) * 100);

                        if (percentProgress > oldPercentProgress)
                        {
                            oldPercentProgress = percentProgress;
                            ReportProgress(new ProgressEventArgs(percentProgress, "Performing classic search... ", NestedIds));
                        }
                    }
                    Interlocked.Add(ref peptideCounter, peptidesInRange);
                };

                if (_singleThreadMode)
                {
                    processProteinRange(0, Proteins.Count);
                }
                else
                {
                    // ClassicSearchEngine's range size: many small ranges, sized from the threads actually used
                    int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;
                    int rangeSize = Math.Max(1, Proteins.Count / (64 * Math.Max(1, maxThreadsPerFile)));
                    Parallel.ForEach(
                        Partitioner.Create(0, Proteins.Count, rangeSize),
                        new ParallelOptions { MaxDegreeOfParallelism = maxThreadsPerFile },
                        range => processProteinRange(range.Item1, range.Item2));
                }
            }

            // only the scans this search changed need resolving
            foreach (int scanIndex in UpdatedIndexes.Keys)
            {
                SpectralMatches[scanIndex]?.ResolveAllAmbiguities();
            }

            HashSet<int> updatedIndexes = UpdatedIndexes.Keys.ToHashSet();
            return new TransientSearchEngineResults(this, updatedIndexes, peptideCounter);
        }

        /// <summary>
        /// The base engine's scan lookup, unless the files' precursor tolerances differ: then the widest window is searched
        /// and narrowed to each scan's file, which needs each interval's centre.
        /// </summary>
        private void FillAcceptableScans(double peptideMonoisotopicMass, List<ScanWithIndexAndNotchInfo> acceptableScans)
        {
            if (!_filterPrecursorPerFile)
            {
                FindAcceptableScans(peptideMonoisotopicMass, acceptableScans);
                return;
            }

            acceptableScans.Clear();
            acceptableScans.AddRange(GetAcceptableScans(peptideMonoisotopicMass, SearchMode));
        }

        private void AddPeptideCandidateToPsm(ScanWithIndexAndNotchInfo scan, double thisScore, IBioPolymerWithSetMods peptide, List<MatchedFragmentIon> matchedIons)
        {
            int scanIndex = scan.ScanIndex;

            if (_singleThreadMode)
            {
                var existingPsm = SpectralMatches[scanIndex];

                if (existingPsm != null)
                {
                    double scoreDiff = thisScore - existingPsm.RunnerUpScore;
                    if (scoreDiff <= -SpectralMatch.ToleranceForScoreDifferentiation)
                        return;
                }

                existingPsm = EnsureWritablePsm(scanIndex, existingPsm);

                UpdatedIndexes.TryAdd(scanIndex, 0);

                if (existingPsm == null)
                {
                    SpectralMatches[scanIndex] = GlobalVariables.AnalyteType == AnalyteType.Oligo
                        ? new OligoSpectralMatch(peptide, scan.Notch, thisScore, scanIndex,
                            ArrayOfSortedMS2Scans[scanIndex], CommonParameters, matchedIons)
                        : new PeptideSpectralMatch(peptide, scan.Notch, thisScore, scanIndex,
                            ArrayOfSortedMS2Scans[scanIndex], CommonParameters, matchedIons);
                }
                else
                {
                    existingPsm.AddOrReplace(peptide, thisScore, scan.Notch,
                        CommonParameters.ReportAllAmbiguity, matchedIons);
                }

                return;
            }

            // Most candidates lose to the match already there, so check without the lock first. The read can
            // be stale, since another thread may be improving this scan's match, but a stale read only means
            // taking the lock for a candidate that then loses the same check below.
            var unlockedPsm = SpectralMatches[scanIndex];
            if (unlockedPsm != null && thisScore - unlockedPsm.RunnerUpScore <= -SpectralMatch.ToleranceForScoreDifferentiation)
                return;

            // The base engine's per-scan locks, as ClassicSearchEngine.AddPeptideCandidateToPsm uses them
            lock (LockFor(scanIndex))
            {
                var existingPsm = SpectralMatches[scanIndex];

                // Double-check after acquiring write lock
                if (existingPsm != null)
                {
                    double scoreDiff = thisScore - existingPsm.RunnerUpScore;
                    if (scoreDiff <= -SpectralMatch.ToleranceForScoreDifferentiation)
                        return;
                }

                existingPsm = EnsureWritablePsm(scanIndex, existingPsm);
                UpdatedIndexes.TryAdd(scanIndex, 0);

                // if the PSM is null, create a new one; otherwise, add or replace the peptide
                if (existingPsm == null)
                {
                    SpectralMatches[scanIndex] = GlobalVariables.AnalyteType == AnalyteType.Oligo
                        ? new OligoSpectralMatch(peptide, scan.Notch, thisScore, scanIndex,
                            ArrayOfSortedMS2Scans[scanIndex], CommonParameters, matchedIons)
                        : new PeptideSpectralMatch(peptide, scan.Notch, thisScore, scanIndex,
                            ArrayOfSortedMS2Scans[scanIndex], CommonParameters, matchedIons);
                }
                else
                {
                    existingPsm.AddOrReplace(peptide, thisScore, scan.Notch,
                        CommonParameters.ReportAllAmbiguity, matchedIons);
                }
            }
        }

        private SpectralMatch EnsureWritablePsm(int scanIndex, SpectralMatch existingPsm)
        {
            if (!_copyOnWriteEnabled || existingPsm == null || UpdatedIndexes.ContainsKey(scanIndex))
            {
                return existingPsm;
            } 
                
            SpectralMatches[scanIndex] = ClonePsmForWrite(existingPsm);
            return SpectralMatches[scanIndex];
        }

        private static SpectralMatch ClonePsmForWrite(SpectralMatch source)
        {
            if (source is PeptideSpectralMatch peptidePsm)
            {
                var bestMatches = peptidePsm.BestMatchingBioPolymersWithSetMods.ToList();
                var cloned = peptidePsm.Clone(bestMatches);
                cloned.PsmFdrInfo = source.PsmFdrInfo?.Clone() ?? new FdrInfo();
                cloned.PeptideFdrInfo = source.PeptideFdrInfo?.Clone();
                return cloned;
            }

            throw new NotSupportedException($"Copy-on-write PSM cloning is not supported for {source.GetType().Name}.");
        }

        private IEnumerable<ScanWithIndexAndNotchInfo> GetAcceptableScans(double peptideMonoisotopicMass, MassDiffAcceptor searchMode)
        {
            // a peptide with an unknown residue (X) has no mass and matches no scan, as in ClassicSearchEngine;
            // a NaN window would otherwise never end the scan loop below and reach every scan
            if (double.IsNaN(peptideMonoisotopicMass))
            {
                yield break;
            }

            foreach (AllowedIntervalWithNotch allowedIntervalWithNotch in searchMode.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(peptideMonoisotopicMass))
            {
                int scanIndex = GetFirstScanWithMassOverOrEqual(allowedIntervalWithNotch.Minimum);
                // the interval is symmetric about the shifted theoretical mass, so its midpoint is that mass
                double centre = (allowedIntervalWithNotch.Minimum + allowedIntervalWithNotch.Maximum) / 2;
                double halfWidthScale = _perFileTolerances?.PrecursorIsPpm == true ? centre : 1.0;
                for (; scanIndex < ArrayOfSortedMS2Scans.Length; scanIndex++)
                {
                    double scanMass = MyScanPrecursorMasses[scanIndex];
                    if (scanMass > allowedIntervalWithNotch.Maximum)
                    {
                        break;
                    }

                    // the interval is the widest of all files; narrow it to this scan's file
                    if (_filterPrecursorPerFile)
                    {
                        if (Math.Abs(scanMass - centre) > halfWidthScale * _scanPrecursorHalfWidth[scanIndex])
                        {
                            continue;
                        }
                    }

                    yield return new ScanWithIndexAndNotchInfo(allowedIntervalWithNotch.Notch, scanIndex);
                }
            }
        }

        private int GetFirstScanWithMassOverOrEqual(double minimum)
        {
            int index = Array.BinarySearch(MyScanPrecursorMasses, minimum);
            if (index < 0)
            {
                index = ~index;
            }

            // index of the first element that is larger than value
            return index;
        }
    }

    public class TransientSearchEngineResults(
        TransientClassicSearchEngine engine,
        HashSet<int> updatedSpectralMatchIndexes,
        int peptidesSearched)
        : MetaMorpheusEngineResults(engine)
    {
        public int PeptidesSearched { get; private set; } = peptidesSearched;
        public HashSet<int> UpdatedSpectralMatchIndexes { get; private set; } = updatedSpectralMatchIndexes;

        public override string ToString()
        {
            var sb = new System.Text.StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("Number of PSMs updated: " + UpdatedSpectralMatchIndexes.Count);
            return sb.ToString();
        }
    }
}
