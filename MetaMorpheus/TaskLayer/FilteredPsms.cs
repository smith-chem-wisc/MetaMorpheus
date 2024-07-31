using Easy.Common.Extensions;
using EngineLayer;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TaskLayer
{
    /// <summary>
    /// Contains a filtered list of PSMs
    /// </summary>
    public class FilteredPsms : IEnumerable<SpectralMatch>
    {
        public List<SpectralMatch> FilteredPsmsList { get; set; }
        /// <summary>
        /// Filter type can have only two values: "q-value" or "pep q-value"
        /// </summary>
        public string FilterType { get; }
        public double FilterThreshold { get; }
        public bool FilteringNotPerformed { get; }
        public bool PeptideLevelFiltering { get; }
        public FilteredPsms(List<SpectralMatch> filteredPsms, string filterType, double filterThreshold, bool filteringNotPerformed, bool peptideLevelFiltering)
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
                case "pep q-value":
                    return psm.GetFdrInfo(PeptideLevelFiltering).PEP_QValue <= FilterThreshold;
                default:
                    return psm.GetFdrInfo(PeptideLevelFiltering).QValue <= FilterThreshold && psm.GetFdrInfo(PeptideLevelFiltering).QValueNotch <= FilterThreshold;
            }
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
            string filterType = "q-value";
            if (pepQValueThreshold < qValueThreshold)
            {
                if (psms.Count() < 100)
                {
                    filteringNotPerformed = true;
                    filterThreshold = 1;
                }
                else
                {
                    filterType = "pep q-value";
                }
            }

            if (!includeHighQValuePsms)
            {
                filteredPsms = filterType.Equals("q-value")
                    ? psms.Where(p => p.GetFdrInfo(filterAtPeptideLevel) != null
                        && p.GetFdrInfo(filterAtPeptideLevel).QValue <= filterThreshold
                        && p.GetFdrInfo(filterAtPeptideLevel).QValueNotch <= filterThreshold).ToList()
                    : psms.Where(p => p.GetFdrInfo(filterAtPeptideLevel) != null && p.GetFdrInfo(filterAtPeptideLevel).PEP_QValue <= filterThreshold).ToList();
            }
            else
            {
                filteredPsms = psms.ToList();
            }

            if (!includeDecoys)
            {
                filteredPsms.RemoveAll(p => p.IsDecoy);
            }
            if (!includeContaminants)
            {
                filteredPsms.RemoveAll(p => p.IsContaminant);
            }
            if (!includeAmbiguous)
            {
                filteredPsms.RemoveAll(p => p.BaseSequence.IsNullOrEmpty());
            }
            if (!includeAmbiguousMods)
            {
                filteredPsms.RemoveAll(p => p.FullSequence.IsNullOrEmpty());
            }
            if (filterAtPeptideLevel)
            {
                //Choose the top scoring PSM for each peptide
                filteredPsms = filteredPsms
                    .OrderByDescending(p => p)
                    .GroupBy(b => b.FullSequence)
                    .Select(b => b.FirstOrDefault()).ToList();
            }

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
}
