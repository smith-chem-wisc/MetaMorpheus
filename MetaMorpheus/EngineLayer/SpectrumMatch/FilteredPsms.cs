using Easy.Common.Extensions;
using EngineLayer.FdrAnalysis;
using System;
using System.Collections.Generic;
using System.Globalization;
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
        /// The filter the matches were judged by: the search's <see cref="CommonParameters.IdentificationFilter"/>, or
        /// the one a caller asked for by name. Under <see cref="FilterType.PepQValue"/>, a match whose PEP was not
        /// trained is judged by its q-value notch instead (<see cref="Passes"/>).
        /// </summary>
        public FilterType FilterType { get; init; }
        public double FilterThreshold { get; init; }
        public bool PeptideLevelFiltering { get; init; }
        /// <summary>
        /// When the PEP q-value filter judged some target matches by their q-value notch, because PEP was not trained
        /// for them, the sentence saying how many and for which digestion agents, for results.txt and the user; null
        /// otherwise.
        /// </summary>
        public string PepFilteringNotDone { get; init; }
        private FilteredPsms(List<SpectralMatch> filteredPsms, FilterType filterType, double filterThreshold, bool peptideLevelFiltering)
        {
            FilteredPsmsList = filteredPsms;
            FilterType = filterType;
            FilterThreshold = filterThreshold;
            PeptideLevelFiltering = peptideLevelFiltering;
        }

        private bool AboveThreshold(SpectralMatch psm) => Passes(psm, PeptideLevelFiltering, FilterType, FilterThreshold);

        /// <summary>
        /// The one identification decision every step of a search reads: whether <paramref name="psm"/> passes
        /// <paramref name="filterType"/> at <paramref name="threshold"/>, at the PSM or the peptide level. A match
        /// without FDR values at that level never passes.
        /// </summary>
        public static bool Passes(SpectralMatch psm, bool peptideLevel, FilterType filterType, double threshold)
        {
            RefuseThePooledQValue(filterType);
            var fdrInfo = psm.GetFdrInfo(peptideLevel);
            return fdrInfo != null && Value(fdrInfo, filterType) <= threshold;
        }

        /// <summary>
        /// The pooled q-value is not an identification filter (<see cref="CommonParameters.IdentificationFilter"/>), so no
        /// caller may ask for it by name either.
        /// </summary>
        private static void RefuseThePooledQValue(FilterType filterType)
        {
            if (filterType is not (FilterType.PepQValue or FilterType.QValueNotch))
                throw new ArgumentException(CommonParameters.NotAFilter(filterType), nameof(filterType));
        }

        /// <summary>
        /// The value a match is judged by. Under the PEP q-value filter, its PEP q-value where PEP was trained for it
        /// (<see cref="FdrInfo.PepTrained"/>), and its q-value notch where it was not. Under the q-value notch filter, its
        /// q-value notch.
        /// </summary>
        private static double Value(FdrInfo fdrInfo, FilterType filterType) => filterType switch
        {
            FilterType.PepQValue => fdrInfo.PepTrained ? fdrInfo.PEP_QValue : fdrInfo.QValueNotch,
            FilterType.QValueNotch => fdrInfo.QValueNotch,
            _ => throw new ArgumentOutOfRangeException(nameof(filterType), filterType, null)
        };

        /// <summary>
        /// The threshold <paramref name="filterType"/> is held to: the PEP q-value threshold for the PEP q-value (and
        /// for the q-value notch of the matches it judges by the notch), the q-value threshold otherwise.
        /// </summary>
        public static double Threshold(CommonParameters commonParams, FilterType filterType) =>
            filterType == FilterType.PepQValue ? commonParams.PepQValueThreshold : commonParams.QValueThreshold;

        /// <summary>The filter and its threshold, as results.txt and the output files state them: "q-value notch &lt;= 0.01".</summary>
        public static string Describe(FilterType filterType, double threshold) =>
            Name(filterType) + " <= " + threshold.ToString(CultureInfo.InvariantCulture);

        /// <summary>The search's own filter and threshold, described as <see cref="Describe(FilterType, double)"/>.</summary>
        public static string Describe(CommonParameters commonParams) =>
            Describe(commonParams.IdentificationFilter, Threshold(commonParams, commonParams.IdentificationFilter));

        private static string Name(FilterType filterType) => filterType switch
        {
            FilterType.PepQValue => "pep q-value",
            FilterType.QValueNotch => "q-value notch",
            FilterType.QValue => "q-value",
            _ => throw new ArgumentOutOfRangeException(nameof(filterType), filterType, null)
        };

        public string GetFilterTypeString() => Name(FilterType);

        /// <summary>
        /// The FDR value this filter decided on for a match, at the level it filtered (as <see cref="AboveThreshold"/>),
        /// e.g. the q-value handed to quantification.
        /// </summary>
        public double GetFilterValue(SpectralMatch psm) => Value(psm.GetFdrInfo(PeptideLevelFiltering), FilterType);

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
        /// Whether <paramref name="psm"/> passes this filter's threshold, at this filter's level (PSM or peptide), by the
        /// rule the PSMs it holds passed. For a match the filter was not given.
        /// </summary>
        public bool PassesThreshold(SpectralMatch psm) => AboveThreshold(psm);

        /// <summary>
        /// Returns the number of PSMs that passed the filtering criteria
        /// </summary>
        public int TargetPsmsAboveThreshold => FilteredPsmsList.Count(psm => !psm.IsDecoy && !psm.IsContaminant && AboveThreshold(psm));

        /// <summary>
        /// Returns a FilteredPsms object that holds every psm that passed the search's identification filter
        /// (<see cref="CommonParameters.IdentificationFilter"/>) at that filter's threshold, or the filter a caller asks
        /// for by name. Under the PEP q-value, a match whose PEP was not trained is judged by its q-value notch at the
        /// same threshold, and <see cref="PepFilteringNotDone"/> says so.
        ///
        /// NOTE: In multiprotease searches, each protease may have its own CommonParameters with different
        /// filtering thresholds. However, this method only accepts a single CommonParameters instance,
        /// so filtering may not be optimal for all proteases when generating a combined filtered PSM list.
        /// </summary>
        /// <param name="psms"> List of spectral match objects to be filtered</param>
        /// <param name="commonParams">The search's parameters: its identification filter and the thresholds</param>
        /// <param name="filterAtPeptideLevel">Filter results at the peptide level (defaults to false) </param>
        /// <param name="includeAmbiguousMods">Keep PSMs whose full sequence did not resolve (two or more full sequences
        /// match the scan)</param>
        /// <param name="includeLocalizationAmbiguousMods">Keep PSMs ambiguous only in where their modifications sit
        /// (<see cref="SpectralMatch.IsLocalizationAmbiguous"/>), whatever <paramref name="includeAmbiguousMods"/> says
        /// about the others. Null, the default, follows <paramref name="includeAmbiguousMods"/>.</param>
        /// <param name="qValueThreshold">Overrides the q-value threshold, which the q-value notch is held to</param>
        /// <param name="pepQValueThreshold">Overrides the PEP q-value threshold, which the PEP q-value is held to</param>
        /// <param name="filterType">The filter, when the caller chooses it. Null (the default) is the search's
        /// <see cref="CommonParameters.IdentificationFilter"/>.</param>
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
            bool? includeLocalizationAmbiguousMods = null,
            FilterType? filterType = null)
        {
            bool includeLocalizationAmbiguous = includeLocalizationAmbiguousMods ?? includeAmbiguousMods;

            // Whether a PSM whose full sequence did not resolve is kept. It is asked whether it is ambiguous only in its
            // modification sites only when the two settings disagree.
            bool KeepUnresolvedFullSequence(SpectralMatch psm) =>
                includeAmbiguousMods == includeLocalizationAmbiguous
                    ? includeAmbiguousMods
                    : psm.IsLocalizationAmbiguous ? includeLocalizationAmbiguous : includeAmbiguousMods;

            FilterType type = filterType ?? commonParams.IdentificationFilter;
            RefuseThePooledQValue(type);
            double threshold = type == FilterType.PepQValue
                ? pepQValueThreshold ?? commonParams.PepQValueThreshold
                : qValueThreshold ?? commonParams.QValueThreshold;

            List<SpectralMatch> candidates = psms.Where(psm =>
                    (includeDecoys || !psm.IsDecoy)
                    && (includeContaminants || !psm.IsContaminant)
                    && (includeAmbiguous || !psm.BaseSequence.IsNullOrEmpty())
                    && (!psm.FullSequence.IsNullOrEmpty() || KeepUnresolvedFullSequence(psm)))
                .ToList();

            List<SpectralMatch> filteredPsms = candidates
                .FilterByQValue(includeHighQValuePsms, threshold, filterAtPeptideLevel, type)
                .CollapseToPeptides(filterAtPeptideLevel)
                .ToList();

            return new FilteredPsms(filteredPsms, type, threshold, filterAtPeptideLevel)
            {
                // At the peptide level, the peptides judged: one match per full sequence
                PepFilteringNotDone = type == FilterType.PepQValue
                    ? PepNotTrainedSentence(candidates.CollapseToPeptides(filterAtPeptideLevel), filterAtPeptideLevel, threshold)
                    : null
            };
        }

        /// <summary>
        /// When PEP was not trained for some target matches at this level, the sentence saying how many, for which
        /// digestion agent (and cleavage-specificity category, which a semi- or non-specific search trains PEP in
        /// separately), and what they were judged by instead. Null when PEP was trained for every one.
        /// </summary>
        private static string PepNotTrainedSentence(IEnumerable<SpectralMatch> candidates, bool peptideLevel, double threshold)
        {
            var untrained = candidates
                .Where(psm => !psm.IsDecoy && psm.GetFdrInfo(peptideLevel) is { PepTrained: false })
                .GroupBy(Group)
                .OrderBy(g => g.Key, StringComparer.Ordinal)
                .Select(g => (g.Key, Count: g.Count()))
                .ToList();
            if (untrained.Count == 0)
            {
                return null;
            }

            int count = untrained.Sum(g => g.Count);
            string plural = count == 1 ? "" : "s";
            string notTrained = peptideLevel
                ? $"peptide-level PEP q-values were not computed for {count} target {GlobalVariables.AnalyteType.GetUniqueFormLabel().ToLower()}{plural}"
                : $"PEP was not trained for {count} target {GlobalVariables.AnalyteType.GetSpectralMatchLabel()}{plural}";
            return $"PEP filtering was requested, but {notTrained} ({string.Join(", ", untrained.Select(g => $"{g.Key}: {g.Count}"))}), " +
                "so those were filtered by " + Describe(FilterType.QValueNotch, threshold) + ".";

            static string Group(SpectralMatch psm)
            {
                var specificity = psm.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer.CleavageSpecificityForFdrCategory;
                string agent = psm.DigestionParams.DigestionAgent.Name;
                return specificity is Omics.Digestion.CleavageSpecificity.Semi or Omics.Digestion.CleavageSpecificity.None
                    ? $"{agent} ({specificity})"
                    : agent;
            }
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

        public static IEnumerable<SpectralMatch> FilterByQValue(this IEnumerable<SpectralMatch> psms, bool includeHighQValuePsms, double qValueThreshold, bool filterAtPeptideLevel, FilterType filterType, bool filterValidated = false)
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
