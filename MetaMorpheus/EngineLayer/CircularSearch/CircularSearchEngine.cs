using EngineLayer.Util;
using MassSpectrometry;
using MzLibUtil;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace EngineLayer.CircularSearch
{
    /// <summary>
    /// A search engine dedicated to head-to-tail cyclic peptides derived from
    /// <see cref="CircularProtein"/> entries.
    ///
    /// DIGESTION
    /// ---------
    /// Each <see cref="CircularProtein"/> is digested via
    /// <see cref="CircularProtein.Digest"/>, which produces:
    ///   • <see cref="CircularPeptideWithSetModifications"/> — intact ring.
    ///     Scored using internal fragment ions only
    ///     (<see cref="CircularPeptideWithSetModifications.FragmentInternally"/>).
    ///   • <see cref="PeptideWithSetModifications"/> — ring-opening linear
    ///     products from a single cut. Scored with BOTH standard terminal ions
    ///     AND internal fragment ions, because these peptides derive from a
    ///     cyclic precursor and internal ions are physically meaningful.
    ///     NOTE: This dual scoring applies ONLY in CircularSearchEngine.
    ///     ClassicSearchEngine scores linear peptides with terminal ions only
    ///     and must not be changed.
    ///
    /// SCORING
    /// -------
    /// Scoring is identical to <see cref="ClassicSearch.ClassicSearchEngine"/>:
    /// matched fragment ions are counted and a PSM score is computed.
    ///
    /// THREAD SAFETY
    /// -------------
    /// One lock object per scan index, matching the ClassicSearchEngine pattern.
    /// </summary>
    public class CircularSearchEngine : MetaMorpheusEngine
    {
        // ── Fields ────────────────────────────────────────────────────────────

        private readonly MassDiffAcceptor SearchMode;
        private readonly List<CircularProtein> CircularProteins;
        private readonly List<Modification> FixedModifications;
        private readonly List<Modification> VariableModifications;
        private readonly SpectralMatch[] SpectralMatches;
        private readonly Ms2ScanWithSpecificMass[] ArrayOfSortedMS2Scans;
        private readonly double[] MyScanPrecursorMasses;
        private readonly int MinInternalFragmentLength;
        private readonly object[] Locks;

        /// <summary>
        /// Internal mutable store for digestion product counts, populated during
        /// <see cref="RunSpecific"/>. External code should use
        /// <see cref="DigestionCountDictionary"/> instead.
        /// </summary>
        private readonly ConcurrentDictionary<(string Accession, string BaseSequence), int>
            _digestionCountDictionary = new();

        /// <summary>
        /// Digestion product counts per (Accession, BaseSequence) pair.
        /// The returned dictionary is a read-only snapshot; it is safe to
        /// enumerate after <see cref="RunSpecific"/> has completed.
        /// </summary>
        public IReadOnlyDictionary<(string Accession, string BaseSequence), int>
            DigestionCountDictionary => _digestionCountDictionary;

        // ── Constructor ───────────────────────────────────────────────────────

        public CircularSearchEngine(
            SpectralMatch[] globalPsms,
            Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans,
            List<Modification> variableModifications,
            List<Modification> fixedModifications,
            List<CircularProtein> circularProteins,
            MassDiffAcceptor searchMode,
            CommonParameters commonParameters,
            List<(string FileName, CommonParameters Parameters)> fileSpecificParameters,
            List<string> nestedIds,
            int minInternalFragmentLength = 3)
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            SpectralMatches = globalPsms;
            ArrayOfSortedMS2Scans = arrayOfSortedMS2Scans;
            MyScanPrecursorMasses = new double[arrayOfSortedMS2Scans.Length];
            for (int i = 0; i < arrayOfSortedMS2Scans.Length; i++)
                MyScanPrecursorMasses[i] = arrayOfSortedMS2Scans[i].PrecursorMass;
            VariableModifications = variableModifications;
            FixedModifications = fixedModifications;
            CircularProteins = circularProteins;
            SearchMode = searchMode;
            MinInternalFragmentLength = minInternalFragmentLength;

            Locks = new object[SpectralMatches.Length];
            for (int i = 0; i < Locks.Length; i++)
                Locks[i] = new object();
        }

        // ── Core search ───────────────────────────────────────────────────────

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Performing circular peptide search...");

            if (!CircularProteins.Any())
                return new MetaMorpheusEngineResults(this);

            int proteinsSearched = 0;
            int oldPercentProgress = 0;

            int maxThreads = CommonParameters.MaxThreadsToUsePerFile;
            var proteinPartitioner = Partitioner.Create(0, CircularProteins.Count);

            Parallel.ForEach(
                proteinPartitioner,
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, _) =>
                {
                    // Per-thread fragment lists, keyed by dissociation type
                    var fragmentsPerDissociationType =
                        new Dictionary<DissociationType, List<Product>>();

                    if (CommonParameters.DissociationType == DissociationType.Autodetect)
                    {
                        foreach (var kv in GlobalVariables.AllSupportedDissociationTypes
                                     .Where(p => p.Value != DissociationType.Autodetect))
                            fragmentsPerDissociationType[kv.Value] = new List<Product>();
                    }
                    else
                    {
                        fragmentsPerDissociationType[CommonParameters.DissociationType] =
                            new List<Product>();
                    }

                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        if (GlobalVariables.StopLoops) return;

                        CircularProtein protein = CircularProteins[i];

                        foreach (var bioPolymer in protein.Digest(
                                     CommonParameters.DigestionParams,
                                     FixedModifications,
                                     VariableModifications))
                        {
                            _digestionCountDictionary.AddOrUpdate(
                                (bioPolymer.Parent.Accession, bioPolymer.BaseSequence),
                                1, (_, v) => v + 1);

                            bool isCircular = bioPolymer is CircularPeptideWithSetModifications;

                            // Determine precursor mass for scan matching.
                            // CircularPeptideWithSetModifications overrides MonoisotopicMass
                            // to subtract water; PeptideWithSetModifications uses base mass.
                            double precursorMass = isCircular
                                ? ((CircularPeptideWithSetModifications)bioPolymer).MonoisotopicMass
                                : ((PeptideWithSetModifications)bioPolymer).MonoisotopicMass;

                            // INVARIANT: Fragment caching contract
                            // ---------------------------------------------------------
                            // All fragment lists are cleared once per bioPolymer so that
                            // stale ions from a previous peptide are never reused.
                            // The computedDissociationTypes set is also reset so each
                            // dissociation type's list is lazily populated exactly once
                            // per bioPolymer on first scan encounter (guarded by the
                            // HashSet.Add check below). This means:
                            //   1. Do NOT move this Clear() inside the scan loop.
                            //   2. Do NOT selectively clear individual keys — every
                            //      key must be reset together.
                            // ---------------------------------------------------------
                            foreach (var fl in fragmentsPerDissociationType.Values)
                                fl.Clear();
                            var computedDissociationTypes = new HashSet<DissociationType>();

                            foreach (var scan in GetAcceptableScans(precursorMass, SearchMode))
                            {
                                var theScan = ArrayOfSortedMS2Scans[scan.ScanIndex];
                                var dissociationType =
                                    CommonParameters.DissociationType == DissociationType.Autodetect
                                        ? theScan.TheScan.DissociationType!.Value
                                        : CommonParameters.DissociationType;

                                if (!fragmentsPerDissociationType.TryGetValue(
                                        dissociationType, out var theorProducts))
                                    continue;

                                // Generate fragments on first encounter for this type
                                if (computedDissociationTypes.Add(dissociationType))
                                {
                                    if (isCircular)
                                    {
                                        // Intact ring: internal fragments only.
                                        // A single backbone cleavage of the intact ring
                                        // produces a linear ion mass-equivalent to the
                                        // precursor — no sequence information is gained
                                        // and these ions are not scored.
                                        ((CircularPeptideWithSetModifications)bioPolymer)
                                            .FragmentInternally(
                                                dissociationType,
                                                MinInternalFragmentLength,
                                                theorProducts);
                                    }
                                    else
                                    {
                                        // Linear ring-opening product: terminal ions AND
                                        // internal ions.
                                        //
                                        // This peptide originated from a cyclic precursor,
                                        // so internal fragment ions are physically meaningful
                                        // and must be scored. This dual scoring is intentional
                                        // and ONLY applies in CircularSearchEngine — do not
                                        // replicate in ClassicSearchEngine.
                                        var linearPeptide = (PeptideWithSetModifications)bioPolymer;

                                        // Terminal ions (b, y, etc.)
                                        linearPeptide.Fragment(
                                            dissociationType,
                                            CommonParameters.DigestionParams.FragmentationTerminus,
                                            theorProducts);

                                        // Internal ions — collected separately because
                                        // FragmentInternally calls Clear() on its output
                                        // list, which would destroy the terminal ions above.
                                        var internalProducts = new List<Product>();
                                        linearPeptide.FragmentInternally(
                                            dissociationType,
                                            MinInternalFragmentLength,
                                            internalProducts);
                                        theorProducts.AddRange(internalProducts);
                                    }
                                }

                                var matchedIons = MatchFragmentIons(
                                    theScan, theorProducts, CommonParameters);

                                double score = CalculatePeptideScore(
                                    theScan.TheScan, matchedIons);

                                AddCandidateToPsm(scan, score, bioPolymer, matchedIons);
                            }
                        }

                        // Progress reporting — thread-safe counter
                        var newCount = Interlocked.Increment(ref proteinsSearched);
                        int pct = (int)((double)newCount / CircularProteins.Count * 100);
                        if (pct > Volatile.Read(ref oldPercentProgress))
                        {
                            // Benign race: two threads may both enter here and
                            // report the same percentage; this is harmless.
                            Volatile.Write(ref oldPercentProgress, pct);
                            ReportProgress(new ProgressEventArgs(
                                pct, "Circular peptide search...", NestedIds));
                        }
                    }
                });

            foreach (var psm in SpectralMatches.Where(p => p != null))
                psm.ResolveAllAmbiguities();

            return new MetaMorpheusEngineResults(this);
        }

        // ── Helpers ───────────────────────────────────────────────────────────

        internal void AddCandidateToPsm(
            ScanWithIndexAndNotchInfo scan,
            double score,
            IBioPolymerWithSetMods peptide,
            List<MatchedFragmentIon> matchedIons)
        {
            if (score < CommonParameters.ScoreCutoff) return;

            lock (Locks[scan.ScanIndex])
            {
                // Compare against RunnerUpScore (not Score) intentionally:
                // candidates scoring between RunnerUpScore and Score must reach
                // AddOrReplace so it can update the runner-up. Candidates slightly
                // below the runner-up are discarded harmlessly inside AddOrReplace.
                // This matches the guard logic in ClassicSearchEngine.
                bool improve = SpectralMatches[scan.ScanIndex] == null
                            || (score - SpectralMatches[scan.ScanIndex].RunnerUpScore)
                               > -SpectralMatch.ToleranceForScoreDifferentiation;

                if (!improve) return;

                if (SpectralMatches[scan.ScanIndex] == null)
                    SpectralMatches[scan.ScanIndex] = new PeptideSpectralMatch(
                        peptide, scan.Notch, score, scan.ScanIndex,
                        ArrayOfSortedMS2Scans[scan.ScanIndex], CommonParameters, matchedIons);
                else
                    SpectralMatches[scan.ScanIndex].AddOrReplace(
                        peptide, score, scan.Notch,
                        CommonParameters.ReportAllAmbiguity, matchedIons);
            }
        }

        private IEnumerable<ScanWithIndexAndNotchInfo> GetAcceptableScans(
            double peptideMass, MassDiffAcceptor searchMode)
        {
            foreach (var interval in searchMode
                         .GetAllowedPrecursorMassIntervalsFromTheoreticalMass(peptideMass))
            {
                int idx = GetFirstScanWithMassOverOrEqual(interval.Minimum);
                while (idx < ArrayOfSortedMS2Scans.Length
                       && MyScanPrecursorMasses[idx] <= interval.Maximum)
                {
                    yield return new ScanWithIndexAndNotchInfo(interval.Notch, idx);
                    idx++;
                }
            }
        }

        private int GetFirstScanWithMassOverOrEqual(double minimum)
        {
            int index = Array.BinarySearch(MyScanPrecursorMasses, minimum);
            return index < 0 ? ~index : index;
        }
    }
}