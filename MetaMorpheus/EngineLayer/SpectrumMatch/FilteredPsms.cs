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
        /// Filter type can have only two values: "q-value" or "pep q-value"
        /// </summary>
        public FilterType FilterType { get; init; }
        public double FilterThreshold { get; init; }
        public bool FilteringNotPerformed { get; init; }
        public bool PeptideLevelFiltering { get; init; }
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
                default:
                    return psm.GetFdrInfo(PeptideLevelFiltering).QValue <= FilterThreshold && psm.GetFdrInfo(PeptideLevelFiltering).QValueNotch <= FilterThreshold;
            }
        }

        public string GetFilterTypeString()
        {
            return FilterType == FilterType.PepQValue ? "pep q-value" : "q-value";
        }

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
        /// Drops each PSM whose full sequence did not resolve (two or more full sequences match the scan) unless
        /// <paramref name="keep"/> accepts it, returning how many were dropped. PSMs with a resolved full sequence are
        /// never dropped. For decisions <see cref="Filter"/> cannot make because they need more than the PSM, such as
        /// its protein groups, or what else was identified in its file.
        /// </summary>
        public int RemoveAmbiguousFullSequences(Func<SpectralMatch, bool> keep)
        {
            return FilteredPsmsList.RemoveAll(psm => psm.FullSequence == null && !keep(psm));
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
        /// <param name="includeAmbiguousMods">Keep PSMs whose full sequence did not resolve (two or more full sequences
        /// match the scan)</param>
        /// <param name="includeLocalizationAmbiguousMods">Keep PSMs ambiguous only in where their modifications sit
        /// (<see cref="SpectralMatch.IsLocalizationAmbiguous"/>), whatever <paramref name="includeAmbiguousMods"/> says
        /// about the others. Null, the default, follows <paramref name="includeAmbiguousMods"/>.</param>
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
            bool? includeLocalizationAmbiguousMods = null)
        {
            bool includeLocalizationAmbiguous = includeLocalizationAmbiguousMods ?? includeAmbiguousMods;

            // Whether a PSM whose full sequence did not resolve is kept. It is asked whether it is ambiguous only in its
            // modification sites only when the two settings disagree.
            bool KeepUnresolvedFullSequence(SpectralMatch psm) =>
                includeAmbiguousMods == includeLocalizationAmbiguous
                    ? includeAmbiguousMods
                    : psm.IsLocalizationAmbiguous ? includeLocalizationAmbiguous : includeAmbiguousMods;

            qValueThreshold ??= commonParams.QValueThreshold;
            pepQValueThreshold ??= commonParams.PepQValueThreshold;
            double filterThreshold = Math.Min((double)qValueThreshold, (double)pepQValueThreshold);
            bool filteringNotPerformed = false;
            List<SpectralMatch> filteredPsms = new List<SpectralMatch>();

            // set the filter type
            FilterType filterType = FilterType.QValue;
            if (pepQValueThreshold < qValueThreshold)
            {
                if (psms.Count() < 100)
                {
                    filteringNotPerformed = true;
                    filterThreshold = 1;
                }
                else
                {
                    filterType = FilterType.PepQValue;
                }
            }

            filteredPsms = psms.Where(psm => 
                    (includeDecoys || !psm.IsDecoy) 
                    && (includeContaminants || !psm.IsContaminant)
                    && (includeAmbiguous || !psm.BaseSequence.IsNullOrEmpty())
                    && (!psm.FullSequence.IsNullOrEmpty() || KeepUnresolvedFullSequence(psm)))
                .FilterByQValue(includeHighQValuePsms, filterThreshold, filterAtPeptideLevel, filterType, true)
                .CollapseToPeptides(filterAtPeptideLevel)
                .ToList();

            return new FilteredPsms(filteredPsms, filterType, filterThreshold, filteringNotPerformed, filterAtPeptideLevel);
        }

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
