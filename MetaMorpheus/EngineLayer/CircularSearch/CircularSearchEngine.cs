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
    ///   • <see cref="CircularPeptideWithSetModifications"/> — intact ring or
    ///     two-cut sub-peptides. Scored using internal fragment ions only
    ///     (<see cref="CircularPeptideWithSetModifications.FragmentInternally"/>).
    ///   • <see cref="PeptideWithSetModifications"/> — ring-opening linear
    ///     products from a single cut. Scored with standard terminal ions.
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
        /// Tracks digestion product counts per protein if enabled.
        /// </summary>
        public readonly ConcurrentDictionary<(string Accession, string BaseSequence), int>
            DigestionCountDictionary = new();

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
            MyScanPrecursorMasses = arrayOfSortedMS2Scans.Select(b => b.PrecursorMass).ToArray();
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

            double proteinsSearched = 0;
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
                            DigestionCountDictionary.AddOrUpdate(
                                (bioPolymer.Parent.Accession, bioPolymer.BaseSequence),
                                1, (_, v) => v + 1);

                            bool isCircular = bioPolymer is CircularPeptideWithSetModifications;

                            // Determine precursor mass for scan matching.
                            // CircularPeptideWithSetModifications overrides MonoisotopicMass
                            // to subtract water; PeptideWithSetModifications uses base mass.
                            double precursorMass = isCircular
                                ? ((CircularPeptideWithSetModifications)bioPolymer).MonoisotopicMass
                                : ((PeptideWithSetModifications)bioPolymer).MonoisotopicMass;

                            // Clear previous fragments
                            foreach (var fl in fragmentsPerDissociationType.Values)
                                fl.Clear();

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
                                if (theorProducts.Count == 0)
                                {
                                    if (isCircular)
                                    {
                                        // Circular: internal fragments only
                                        ((CircularPeptideWithSetModifications)bioPolymer)
                                            .FragmentInternally(
                                                dissociationType,
                                                MinInternalFragmentLength,
                                                theorProducts);
                                    }
                                    else
                                    {
                                        // Linear ring-opening: standard terminal ions
                                        bioPolymer.Fragment(
                                            dissociationType,
                                            CommonParameters.DigestionParams.FragmentationTerminus,
                                            theorProducts);
                                    }
                                }

                                var matchedIons = MatchFragmentIons(
                                    theScan, theorProducts, CommonParameters);

                                double score = CalculatePeptideScore(
                                    theScan.TheScan, matchedIons);

                                AddCandidateToPsm(scan, score, bioPolymer, matchedIons);
                            }
                        }

                        // Progress reporting
                        proteinsSearched++;
                        int pct = (int)(proteinsSearched / CircularProteins.Count * 100);
                        if (pct > oldPercentProgress)
                        {
                            oldPercentProgress = pct;
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