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
    /// Contains a filtered list of PSMs. Generated within the PostSearchAnalysisTask
    /// </summary>
    public class FilteredPsms : IEnumerable<SpectralMatch>
    {
        public List<SpectralMatch> Psms { get; }
        /// <summary>
        /// Filter type can have only two values: "q-value" or "pep q-value"
        /// </summary>
        public string FilterType { get; }
        public double FilterThreshold { get; }
        public bool FilteringNotPerformed { get; }
        public bool PeptideLevelFiltering { get; }
        public FilteredPsms(List<SpectralMatch> psms, string filterType, double filterThreshold, bool filteringNotPerformed, bool peptideLevelFiltering)
        {
            Psms = psms;
            FilterType = filterType;
            FilterThreshold = filterThreshold;
            FilteringNotPerformed = filteringNotPerformed;
            PeptideLevelFiltering = peptideLevelFiltering;
        }

        private bool AboveThreshold(SpectralMatch psm)
        {
            switch(FilterType)
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
        public int PsmsAboveThreshold => Psms.Count(psm => AboveThreshold(psm));

        public IEnumerator<SpectralMatch> GetEnumerator()
        {
            return Psms.GetEnumerator();
        }

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
            return Psms.GetEnumerator();
        }
    }

    public class PsmFilter
    {
        CommonParameters CommonParams { get; }

        public PsmFilter(CommonParameters commonParameters)
        {
            CommonParams = commonParameters;
        }


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
        public FilteredPsms Filter(List<SpectralMatch> psms,
            bool includeDecoys = true,
            bool includeContaminants = true,
            bool includeAmbiguous = false,
            bool includeAmbiguousMods = true,
            bool includeHighQValuePsms = false,
            double? qValueThreshold = null,
            double? pepQValueThreshold = null,
            bool filterAtPeptideLevel = false)
        {

            qValueThreshold ??= CommonParams.QValueThreshold;
            pepQValueThreshold ??= CommonParams.PepQValueThreshold;
            double filterThreshold = Math.Min((double)qValueThreshold, (double)pepQValueThreshold);
            bool filteringNotPerformed = false;
            List<SpectralMatch> filteredPsms = new List<SpectralMatch>();

            // set the filter type
            string filterType = "q-value";
            if (pepQValueThreshold < qValueThreshold)
            {
                if (psms.Count < 100)
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
                    ? psms.Where(p => p.GetFdrInfo(filterAtPeptideLevel).QValue <= filterThreshold && p.GetFdrInfo(filterAtPeptideLevel).QValueNotch <= filterThreshold).ToList()
                    : psms.Where(p => p.GetFdrInfo(filterAtPeptideLevel).PEP_QValue <= filterThreshold)
                    .ToList();
            }
            else
            {
                filteredPsms = psms;
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
                filteredPsms = filteredPsms
                    .GroupBy(b => b.FullSequence)
                    .Select(b => b.FirstOrDefault()).ToList();
            }

            return new FilteredPsms(filteredPsms, filterType, filterThreshold, filteringNotPerformed, filterAtPeptideLevel);
        }
    }
}
