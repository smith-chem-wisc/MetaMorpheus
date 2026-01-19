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
        public FilteredPsms(List<SpectralMatch> filteredPsms, FilterType filterType, double filterThreshold, bool filteringNotPerformed, bool peptideLevelFiltering)
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
        /// Returns the number of PSMs that passed the filtering criteria
        /// </summary>
        public int TargetPsmsAboveThreshold => FilteredPsmsList.Count(psm => !psm.IsDecoy && !psm.IsContaminant && AboveThreshold(psm));

        /// <summary>
        /// Returns a FilteredPsms object that holds every psm that passed the filtering criteria.
        /// Q-Value and PEP Q-Value thresholds are read from common parameters by default, but can be overridden
        /// Q-Value and PEP Q-Value filtering are mutually exculsive.
        /// In cases where PEP filtering was selected but PEP wasn't performed due to insufficient PSMs, 
        /// filtering defaults to Q and Q_Notch.
        /// </summary>
        /// <param name="psms"> List of spectral match objects to be filtered</param>
        /// <param name="filterAtPeptideLevel">Filter results at the peptide level (defaults to false) </param>
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
            bool filterAtPeptideLevel = false)
        {

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
                    && (includeAmbiguousMods || !psm.FullSequence.IsNullOrEmpty()))
                .FilterByQValue(includeHighQValuePsms, filterThreshold, filterAtPeptideLevel, filterType)
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

    public static class FilteredPsmsExtensions
    {
        public static IEnumerable<SpectralMatch> CollapseToPeptides(this IEnumerable<SpectralMatch> psms, bool filterAtPeptideLevel)
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

        public static IEnumerable<SpectralMatch> FilterByQValue(this IEnumerable<SpectralMatch> psms, bool includeHighQValuePsms, double qValueThreshold, bool filterAtPeptideLevel, FilterType filterType)
        {
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
