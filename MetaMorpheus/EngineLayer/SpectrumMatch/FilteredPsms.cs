using Easy.Common.Extensions;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.SpectrumMatch
{
    /// <summary>
    /// Contains a filtered list of PSMs.
    /// All properties within this class are read-only, and should only be set on object construction
    /// </summary>
    public class FilteredPsms : IEnumerable<SpectralMatch>
    {
        public List<SpectralMatch> FilteredPsmsList { get; private set; }
        /// <summary>
        /// The value the matches were filtered on: "q-value", "pep q-value" or, only when a caller asked for it,
        /// "q-value notch"
        /// </summary>
        public FilterType FilterType { get; init; }
        public double FilterThreshold { get; init; }
        public bool FilteringNotPerformed { get; init; }
        public bool PeptideLevelFiltering { get; init; }
        /// <summary>
        /// When PEP filtering was requested (a PEP q-value threshold below the q-value threshold) but not done, the
        /// sentence saying what was done instead, for results.txt and the user; null otherwise.
        /// </summary>
        public string PepFilteringNotDone { get; init; }
        private FilteredPsms(List<SpectralMatch> filteredPsms, FilterType filterType, double filterThreshold, bool filteringNotPerformed, bool peptideLevelFiltering)
        {
            FilteredPsmsList = filteredPsms;
            FilterType = filterType;
            FilterThreshold = filterThreshold;
            FilteringNotPerformed = filteringNotPerformed;
            PeptideLevelFiltering = peptideLevelFiltering;
        }

        private bool AboveThreshold(SpectralMatch psm)
        {
            if (psm.GetFdrInfo(PeptideLevelFiltering) == null) return false;

            switch (FilterType)
            {
                case FilterType.PepQValue:
                    return psm.GetFdrInfo(PeptideLevelFiltering).PEP_QValue <= FilterThreshold;
                case FilterType.QValueNotch:
                    return psm.GetFdrInfo(PeptideLevelFiltering).QValueNotch <= FilterThreshold;
                default:
                    return psm.GetFdrInfo(PeptideLevelFiltering).QValue <= FilterThreshold && psm.GetFdrInfo(PeptideLevelFiltering).QValueNotch <= FilterThreshold;
            }
        }

        public string GetFilterTypeString() => FilterType switch
        {
            FilterType.PepQValue => "pep q-value",
            FilterType.QValueNotch => "q-value notch",
            _ => "q-value"
        };

        /// <summary>
        /// The FDR value this filter decided on for a match, at the level it filtered (as <see cref="AboveThreshold"/>),
        /// e.g. the q-value handed to quantification.
        /// </summary>
        public double GetFilterValue(SpectralMatch psm) => FilterType switch
        {
            FilterType.PepQValue => psm.GetFdrInfo(PeptideLevelFiltering).PEP_QValue,
            FilterType.QValueNotch => psm.GetFdrInfo(PeptideLevelFiltering).QValueNotch,
            _ => psm.GetFdrInfo(PeptideLevelFiltering).QValue
        };

        /// <summary>
        /// This method should only be called when filtered PSMs are modified for the purpose of SILAC analysis
        /// </summary>
        /// <param name="silacPsms"></param>
        public void SetSilacFilteredPsms(List<SpectralMatch> silacPsms)
        {
            FilteredPsmsList = silacPsms;
        }

        /// <summary>
        /// Drops PSMs whose monoisotopic mass did not resolve, returning how many were dropped.
        /// A resolved sequence does not imply a resolved mass, so callers that require a mass
        /// (quantification) must screen for it separately.
        /// </summary>
        public int RemovePsmsWithoutResolvedMass()
        {
            return FilteredPsmsList.RemoveAll(psm => !psm.BioPolymerWithSetModsMonoisotopicMass.HasValue);
        }

        /// <summary>
        /// Returns the number of PSMs that passed the filtering criteria
        /// </summary>
        public int TargetPsmsAboveThreshold => FilteredPsmsList.Count(psm => !psm.IsDecoy && !psm.IsContaminant && AboveThreshold(psm));

        /// <summary>
        /// Returns a FilteredPsms object that holds every psm that passed the filtering criteria.
        /// Q-Value and PEP Q-Value thresholds are read from common parameters by default, but can be overridden
        /// Q-Value and PEP Q-Value filtering are mutually exclusive.
        /// In cases where PEP filtering was selected but PEP wasn't performed due to insufficient PSMs, 
        /// filtering defaults to Q and Q_Notch.
        /// 
        /// NOTE: In multiprotease searches, each protease may have its own CommonParameters with different 
        /// filtering thresholds. However, this method only accepts a single CommonParameters instance, 
        /// so filtering may not be optimal for all proteases when generating a combined filtered PSM list.
        /// Consider passing explicit qValueThreshold and pepQValueThreshold parameters to ensure consistent 
        /// filtering across all proteases.
        /// </summary>
        /// <param name="psms"> List of spectral match objects to be filtered</param>
        /// <param name="commonParams">Common parameters containing default filter thresholds. In multiprotease scenarios, 
        /// consider using explicit threshold parameters instead.</param>
        /// <param name="filterAtPeptideLevel">Filter results at the peptide level (defaults to false) </param>
        /// <param name="filterType">The value to filter on, chosen by the caller. Null (the default) infers it from the
        /// thresholds, as MetaMorpheus always has. Given, it is used as asked: <see cref="FilterType.PepQValue"/> at the
        /// PEP q-value threshold, <see cref="FilterType.QValue"/> and <see cref="FilterType.QValueNotch"/> at the q-value
        /// threshold, with no count rule and no fallback. A PEP q-value filter on matches that carry no PEP q-value
        /// throws instead of returning nothing.</param>
        /// <returns> A FilteredPsms object</returns>
        public static FilteredPsms Filter(IEnumerable<SpectralMatch> psms,
            CommonParameters commonParams,
            bool includeDecoys = true,
            bool includeContaminants = true,
            bool includeAmbiguous = false,
            bool includeAmbiguousMods = true,
            bool includeHighQValuePsms = false,
            double? qValueThreshold = null,
            double? pepQValueThreshold = null,
            bool filterAtPeptideLevel = false,
            FilterType? filterType = null)
        {

            qValueThreshold ??= commonParams.QValueThreshold;
            pepQValueThreshold ??= commonParams.PepQValueThreshold;
            double filterThreshold = Math.Min((double)qValueThreshold, (double)pepQValueThreshold);
            bool filteringNotPerformed = false;
            string pepFilteringNotDone = null;

            List<SpectralMatch> candidates = psms.Where(psm =>
                    (includeDecoys || !psm.IsDecoy)
                    && (includeContaminants || !psm.IsContaminant)
                    && (includeAmbiguous || !psm.BaseSequence.IsNullOrEmpty())
                    && (includeAmbiguousMods || !psm.FullSequence.IsNullOrEmpty()))
                .ToList();
            string label = GlobalVariables.AnalyteType.GetSpectralMatchLabel();

            // set the filter type
            FilterType type = FilterType.QValue;
            if (filterType is { } requested)
            {
                type = requested;
                filterThreshold = requested == FilterType.PepQValue ? (double)pepQValueThreshold : (double)qValueThreshold;
                if (requested == FilterType.PepQValue && !includeHighQValuePsms && candidates.Count > 0 && NoPepQValues(candidates, filterAtPeptideLevel))
                    throw new MetaMorpheusException($"PEP q-value filtering was requested, but none of these {candidates.Count} {label}s has a PEP q-value: PEP was not computed for them.");
            }
            else if (pepQValueThreshold < qValueThreshold)
            {
                int inputCount = psms.Count();
                if (inputCount < 100)
                {
                    filteringNotPerformed = true;
                    filterThreshold = 1;
                    // The count rule counts the input, decoys included, as it always has.
                    pepFilteringNotDone = $"PEP filtering was requested, but there were fewer than 100 {label}s ({inputCount}, counting decoys), so results were not filtered.";
                }
                else if (candidates.Count > 0 && NoPepQValues(candidates, filterAtPeptideLevel))
                {
                    // PEP is trained only when a digestion agent has more than 1000 matches (FdrAnalysisEngine), so every
                    // PEP q-value here is still the initial 2. A PEP q-value filter would drop every match. At the peptide
                    // level this also happens with PSM-level PEP trained: more than 1000 PSMs but at most 1000 peptides.
                    pepFilteringNotDone = filterAtPeptideLevel
                        ? "PEP filtering was requested, but PEP q-values were not computed at the peptide level, so results were filtered by q-value."
                        : $"PEP filtering was requested, but PEP was not computed for these {label}s, so results were filtered by q-value.";
                }
                else
                {
                    type = FilterType.PepQValue;
                }
            }

            List<SpectralMatch> filteredPsms = candidates
                .FilterByQValue(includeHighQValuePsms, filterThreshold, filterAtPeptideLevel, type, true)
                .CollapseToPeptides(filterAtPeptideLevel)
                .ToList();

            return new FilteredPsms(filteredPsms, type, filterThreshold, filteringNotPerformed, filterAtPeptideLevel)
            {
                PepFilteringNotDone = pepFilteringNotDone
            };
        }

        /// <summary>
        /// True when no match carries a computed PEP q-value at this level: each is still FdrInfo's initial 2.
        /// </summary>
        private static bool NoPepQValues(List<SpectralMatch> psms, bool filterAtPeptideLevel) =>
            psms.All(p => p.GetFdrInfo(filterAtPeptideLevel) != null && Math.Abs(p.GetFdrInfo(filterAtPeptideLevel).PEP_QValue - 2) < 0.1);

        public IEnumerator<SpectralMatch> GetEnumerator()
        {
            return FilteredPsmsList.GetEnumerator();
        }

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
            return FilteredPsmsList.GetEnumerator();
        }
    }

    internal static class FilteredPsmsExtensions
    {
        internal static IEnumerable<SpectralMatch> CollapseToPeptides(this IEnumerable<SpectralMatch> psms, bool filterAtPeptideLevel)
        {
            if (!filterAtPeptideLevel)
            {
                return psms;
            }
            else
            {
                return psms
                    .OrderBy(p => p.FdrInfo.PEP) // Order by PEP first
                    .ThenByDescending(p => p) // Use the default comparer to break ties
                    .GroupBy(p => p.FullSequence) 
                    .Select(p => p.FirstOrDefault()) // Choose the PSM with the lowest PEP for each peptide
                    .OrderByDescending(p => p); // Use the default comparer to order the resulting list, so that the output is ordered by MetaMorpheus Score
            }
        }

        internal static IEnumerable<SpectralMatch> FilterByQValue(this IEnumerable<SpectralMatch> psms, bool includeHighQValuePsms, double qValueThreshold, bool filterAtPeptideLevel, FilterType filterType, bool filterValidated = false)
        {
            // If we validate the filter ahead of time, no need for a second iteration through the PSMs to check for the presence of PEP-QValues values. 
            if (!filterValidated &&  filterType == FilterType.PepQValue && psms.All(p => p.GetFdrInfo(filterAtPeptideLevel) != null && Math.Abs(p.GetFdrInfo(filterAtPeptideLevel).PEP_QValue - 2) < 0.1 ))
            {
                filterType = FilterType.QValue;
            }

            foreach (var psm in psms)
            {
                if (includeHighQValuePsms)
                {
                    yield return psm;
                }
                else if (psm.GetFdrInfo(filterAtPeptideLevel) != null)
                {
                    switch(filterType)
                    {
                        case FilterType.PepQValue:
                            if (psm.GetFdrInfo(filterAtPeptideLevel).PEP_QValue <= qValueThreshold)
                            {
                                yield return psm;
                            }
                            break;
                        case FilterType.QValueNotch:
                            if (psm.GetFdrInfo(filterAtPeptideLevel).QValueNotch <= qValueThreshold)
                            {
                                yield return psm;
                            }
                            break;
                        case FilterType.QValue:
                        default:
                            if (psm.GetFdrInfo(filterAtPeptideLevel).QValue <= qValueThreshold && psm.GetFdrInfo(filterAtPeptideLevel).QValueNotch <= qValueThreshold)
                            {
                                yield return psm;
                            }
                            break;
                    }
                }
            }
        }
    }
}
