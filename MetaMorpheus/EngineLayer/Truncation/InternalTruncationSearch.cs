using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MassSpectrometry;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer.Truncation
{
    /// <summary>
    /// Direct internal-fragment search over the parent-restricted proteoform set — the design ProSightPD's
    /// subsequence search uses, restricted to confident parents. For each scan it enumerates internal spans of
    /// each parent whose mass matches the precursor — found by prefix-sum arithmetic + binary search so only the
    /// few mass-matches are actually built — scores each candidate's OWN b/y ladder with the
    /// fragmentation-propensity-weighted score (<see cref="FragmentationPropensity"/>, the C-Score core), and
    /// keeps the best candidate that clears a strong bilateral gate (≥ minIonsPerTerminus matched ions on EACH
    /// new terminus). Internal PSMs go to a SEPARATE (per-class) FDR — they are a different score regime from
    /// terminal/intact matches, and pooling them against the easy intact threshold is what sank the prior attempt.
    ///
    /// This is necessarily distinct from Pass 2: an internal fragment shares NEITHER ion ladder with its
    /// full-length parent, so the parent fragment index cannot find it; here the candidate's own fragments are
    /// scored. Parallelized across scans (the parent prefix sums are precomputed once, also in parallel).
    /// </summary>
    public static class InternalTruncationSearch
    {
        private const double WaterMonoisotopicMass = 18.0105646;
        // Prefix-window half-width for the internal-span mass search, sized to cover up to ±3 monoisotopic
        // notches plus precursor tolerance. NOTE: this is a fixed constant, NOT derived from the supplied
        // acceptor — if an acceptor permitting a wider notch span than this is ever used, some valid internal
        // spans could be missed. The shipped acceptors stay within ±3 notches; widen this constant (or derive
        // it from the acceptor's notch range) before pairing the internal search with a wider acceptor. (#8)
        private const double NotchSearchWindow = 3.6;

        private static readonly double[] ResidueMass = BuildResidueTable();

        private static double[] BuildResidueTable()
        {
            var t = new double[128];
            for (int i = 0; i < t.Length; i++) t[i] = double.NaN;
            t['G'] = 57.02146; t['A'] = 71.03711; t['S'] = 87.03203; t['P'] = 97.05276; t['V'] = 99.06841;
            t['T'] = 101.04768; t['C'] = 103.00919; t['L'] = 113.08406; t['I'] = 113.08406; t['N'] = 114.04293;
            t['D'] = 115.02694; t['Q'] = 128.05858; t['K'] = 128.09496; t['E'] = 129.04259; t['M'] = 131.04049;
            t['H'] = 137.05891; t['F'] = 147.06841; t['R'] = 156.10111; t['Y'] = 163.06333; t['W'] = 186.07931;
            return t;
        }

        /// <summary>Prefix sums of residue+residue-mod masses for one parent (Prefix[k] = mass of residues 0..k-1);
        /// null if the parent contains an unknown residue (skipped).</summary>
        private static double[] BuildPrefix(TruncationParent parent)
        {
            string sequence = parent.Proteoform.BaseSequence;
            int length = sequence.Length;
            var prefix = new double[length + 1];
            IReadOnlyDictionary<int, Modification> mods = parent.Proteoform.AllModsOneIsNterminus;
            for (int i = 0; i < length; i++)
            {
                char residue = sequence[i];
                double mass = residue < 128 ? ResidueMass[residue] : double.NaN;
                if (double.IsNaN(mass)) return null;
                if (mods != null && mods.TryGetValue(i + 2, out Modification mod) && mod?.MonoisotopicMass != null)
                {
                    mass += mod.MonoisotopicMass.Value; // residue i (0-based) -> mod key i+2
                }
                prefix[i + 1] = prefix[i] + mass;
            }
            return prefix;
        }

        /// <summary>
        /// Enumerates internal fragments directly: for every scan, every parent, and every start residue, the
        /// prefix array is binary-searched for the ends whose span mass lands in the notch window; each
        /// surviving span is built, fragmented and matched.
        /// </summary>
        /// <remarks>
        /// COST. There is no fragment index here, so the work is proportional to (parent count x parent
        /// length) per scan -- it is the PARENT SET that sets the price, not the database the parents came
        /// from. Measured on this branch rather than asserted (Release, one scan so the outer Parallel.For
        /// over scans is not hiding anything, synthetic repeat-alphabet parents, exact-notch acceptor):
        ///
        ///   * linear in parent count -- 0.023 ms per 50-residue parent, flat from 25 through 800 parents;
        ///   * linear in parent length -- 0.111, 0.192 and 0.386 ms per parent at 100, 200 and 400 residues.
        ///
        /// That is about 1 microsecond per residue, per parent, per scan on one thread. Extrapolating to
        /// 5,000 MS2 scans against parents averaging 350 residues, on 60 threads: roughly 350 ms per scan and
        /// half a minute in total against a tag-filtered set of ~1,000 parents, and roughly 7 s per scan and
        /// ten minutes in total against an unfiltered ~20,000-parent human set.
        ///
        /// So the sequence-tag filter buys about a factor of twenty, but this search is not GATED on it --
        /// the unfiltered case is minutes, not hours, on a many-core box. The default stays off for the FDR
        /// reason documented on the parameter, not for cost. Synthetic parents make these a floor: real
        /// proteoforms carry modifications, which multiply the forms per parent.
        /// </remarks>
        public static List<TruncationPsm> Run(IReadOnlyList<TruncationParent> parents, Ms2ScanWithSpecificMass[] scans,
            CommonParameters commonParameters, MassDiffAcceptor chopAcceptor, int minIonsPerTerminus, double maxParentMass)
        {
            // Precompute prefix sums per parent once (parallel); skip oversized / unknown-residue parents.
            var prefixes = new double[parents.Count][];
            Parallel.For(0, parents.Count, i =>
            {
                if (parents[i].MonoisotopicMass <= maxParentMass) prefixes[i] = BuildPrefix(parents[i]);
            });

            var results = new TruncationPsm[scans.Length];
            var parallelOptions = new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, commonParameters.MaxThreadsToUsePerFile) };

            Parallel.For(0, scans.Length, parallelOptions, scanIndex =>
            {
                Ms2ScanWithSpecificMass scan = scans[scanIndex];
                double precursorMass = scan.PrecursorMass;
                double totalIntensity = scan.TheScan.TotalIonCurrent;
                DissociationType dissociationType = TruncationSearchEngine.ResolveDissociationType(commonParameters, scan);
                var products = new List<Product>();
                TruncationPsm best = null;

                for (int pi = 0; pi < parents.Count; pi++)
                {
                    double[] prefix = prefixes[pi];
                    if (prefix == null) continue;
                    TruncationParent parent = parents[pi];
                    if (parent.MonoisotopicMass < precursorMass) continue; // an internal fragment is lighter than its parent
                    int length = prefix.Length - 1;

                    for (int start = 1; start <= length - 2; start++)
                    {
                        // span [start..end] mass = water + prefix[end+1] - prefix[start]; want it within the notch window of M
                        double lowPrefix = (precursorMass - NotchSearchWindow) - WaterMonoisotopicMass + prefix[start];
                        double highPrefix = (precursorMass + NotchSearchWindow) - WaterMonoisotopicMass + prefix[start];
                        if (prefix[length - 1] < lowPrefix) break; // even the largest end is too light here; larger starts only worse

                        for (int k = LowerBound(prefix, start + 1, length, lowPrefix); k <= length - 1 && prefix[k] <= highPrefix; k++)
                        {
                            int end = k - 1;
                            PeptideWithSetModifications form = ProteoformChopper.BuildTruncatedForm(
                                parent.Proteoform, chopFromN: start, chopFromC: (length - 1) - end);

                            int notch = chopAcceptor.Accepts(precursorMass, form.MonoisotopicMass);
                            if (notch < 0) continue;

                            form.Fragment(dissociationType, FragmentationTerminus.Both, products, commonParameters.FragmentationParameters);
                            List<MatchedFragmentIon> matched = MetaMorpheusEngine.MatchFragmentIons(scan, products, commonParameters);

                            int nIons = 0, cIons = 0;
                            foreach (MatchedFragmentIon ion in matched)
                            {
                                if (ion.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.N) nIons++;
                                else if (ion.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.C) cIons++;
                            }
                            if (nIons < minIonsPerTerminus || cIons < minIonsPerTerminus) continue; // strong bilateral gate

                            double score = FragmentationPropensity.Score(matched, form.BaseSequence, totalIntensity, dissociationType);
                            if (best == null || score > best.Score)
                            {
                                var psm = new PeptideSpectralMatch(form, notch, score, scanIndex, scan, commonParameters, matched);
                                best = new TruncationPsm
                                {
                                    SpectralMatch = psm,
                                    TruncatedForm = form,
                                    TruncationProductType = TruncationPass3.InternalTruncation,
                                    ProteinAccessions = parent.ProteinAccession,
                                    ScanIndex = scanIndex
                                };
                            }
                        }
                    }
                }

                results[scanIndex] = best;
            });

            return results.Where(r => r != null).ToList();
        }

        /// <summary>First index in [lo, hi) whose prefix value is ≥ <paramref name="value"/> (prefix is ascending).</summary>
        private static int LowerBound(double[] prefix, int lo, int hi, double value)
        {
            while (lo < hi)
            {
                int mid = (lo + hi) >> 1;
                if (prefix[mid] < value) lo = mid + 1;
                else hi = mid;
            }
            return lo;
        }
    }
}
