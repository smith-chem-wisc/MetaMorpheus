using Easy.Common.Extensions;
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
        /// <summary>
        /// The identification filter the search named (<see cref="CommonParameters.IdentificationFilter"/>), or
        /// <see cref="IdentificationFilterMode.Thresholds"/> when the filter was inferred from the thresholds or the
        /// caller named a <see cref="SpectrumMatch.FilterType"/>. When named, <see cref="IdentificationRule"/> judges
        /// every match, strictly below <see cref="FilterThreshold"/>, and <see cref="FilterType"/> is the value the
        /// matches were judged on (the notch when digestion agents differ).
        /// </summary>
        public IdentificationFilterMode IdentificationFilter { get; init; }
        /// <summary>
        /// Under a named identification filter, the values the matches were judged on, e.g. "pep q-value", or
        /// "pep q-value or q-value notch" when digestion agents differ; null otherwise.
        /// </summary>
        private string JudgedOn { get; init; }
        private bool IsExplicit => IdentificationFilter != IdentificationFilterMode.Thresholds;
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
            if (IsExplicit) return IdentificationRule.Passes(psm, PeptideLevelFiltering, IdentificationFilter, FilterThreshold);

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

        /// <summary>True when the match passes this filter, e.g. the mzIdentML passThreshold.</summary>
        public bool Passes(SpectralMatch psm) => AboveThreshold(psm);

        public string GetFilterTypeString() => JudgedOn ?? GetFilterTypeString(FilterType);

        public static string GetFilterTypeString(FilterType filterType) => filterType switch
        {
            FilterType.PepQValue => "pep q-value",
            FilterType.QValueNotch => "q-value notch",
            FilterType.QValue => "q-value",
            _ => throw new ArgumentOutOfRangeException(nameof(filterType), filterType, null)
        };

        /// <summary>
        /// The comparison and threshold for summary text: "&lt;= 0.01" when the filter was inferred from the
        /// thresholds (unchanged), "&lt; 0.01" under a named identification filter, whose comparison is strict.
        /// </summary>
        public string GetThresholdString() => IsExplicit
            ? "< " + FilterThreshold.ToString(CultureInfo.InvariantCulture)
            : "<= " + Math.Round(FilterThreshold, 2);

        /// <summary>
        /// Applies this filter's comparison to a value already reduced to a q-value, e.g. a protein group's: strict
        /// under a named identification filter, inclusive otherwise.
        /// </summary>
        public bool PassesThreshold(double value) => IsExplicit ? value < FilterThreshold : value <= FilterThreshold;

        /// <summary>
        /// The FDR value this filter decided on for a match, at the level it filtered (as <see cref="AboveThreshold"/>),
        /// e.g. the q-value handed to quantification.
        /// </summary>
        public double GetFilterValue(SpectralMatch psm) => IsExplicit
            ? IdentificationRule.GetValue(psm, PeptideLevelFiltering, IdentificationFilter)
            : FilterType switch
            {
                FilterType.PepQValue => psm.GetFdrInfo(PeptideLevelFiltering).PEP_QValue,
                FilterType.QValueNotch => psm.GetFdrInfo(PeptideLevelFiltering).QValueNotch,
                FilterType.QValue => psm.GetFdrInfo(PeptideLevelFiltering).QValue,
                _ => throw new ArgumentOutOfRangeException(nameof(FilterType), FilterType, null)
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
        /// <param name="filterType">The value to filter on, chosen by the caller. Null (the default) infers it from the
        /// thresholds, as MetaMorpheus always has. Given, it is used as asked: <see cref="FilterType.PepQValue"/> at the
        /// PEP q-value threshold, <see cref="FilterType.QValue"/> and <see cref="FilterType.QValueNotch"/> at the q-value
        /// threshold, with no count rule. A PEP q-value filter on matches that carry no PEP q-value filters by q-value
        /// at the PEP q-value threshold instead of dropping every match, and says so in
        /// <see cref="PepFilteringNotDone"/>, so a search still writes its results.</param>
        /// <remarks>When <paramref name="commonParams"/> names an identification filter
        /// (<see cref="CommonParameters.IdentificationFilter"/>) and no <paramref name="filterType"/> is given, that filter
        /// decides, through <see cref="IdentificationRule"/>, at its threshold (computed from the two thresholds above).
        /// Under <see cref="IdentificationFilterMode.PepQValue"/>, matches of a digestion agent PEP was not trained for
        /// throw, unless <paramref name="includeHighQValuePsms"/> keeps every match anyway.</remarks>
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

            qValueThreshold ??= commonParams.QValueThreshold;
            pepQValueThreshold ??= commonParams.PepQValueThreshold;
            double filterThreshold = Math.Min((double)qValueThreshold, (double)pepQValueThreshold);
            bool filteringNotPerformed = false;
            string pepFilteringNotDone = null;

            List<SpectralMatch> candidates = psms.Where(psm =>
                    (includeDecoys || !psm.IsDecoy)
                    && (includeContaminants || !psm.IsContaminant)
                    && (includeAmbiguous || !psm.BaseSequence.IsNullOrEmpty())
                    && (!psm.FullSequence.IsNullOrEmpty() || KeepUnresolvedFullSequence(psm)))
                .ToList();
            string label = GlobalVariables.AnalyteType.GetSpectralMatchLabel();

            if (filterType == null && IdentificationRule.IsExplicit(commonParams))
            {
                return FilterByName(candidates, commonParams, includeHighQValuePsms, (double)qValueThreshold, (double)pepQValueThreshold, filterAtPeptideLevel);
            }

            // set the filter type
            FilterType type = FilterType.QValue;
            if (filterType is { } requested)
            {
                type = requested;
                filterThreshold = requested == FilterType.PepQValue ? (double)pepQValueThreshold : (double)qValueThreshold;
                if (requested == FilterType.PepQValue && candidates.Count > 0 && NoPepQValues(candidates, filterAtPeptideLevel))
                {
                    type = FilterType.QValue;
                    pepFilteringNotDone = PepNotComputedSentence(filterAtPeptideLevel, label);
                }
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
                    pepFilteringNotDone = PepNotComputedSentence(filterAtPeptideLevel, label);
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
        /// The filter the search named: each match judged by <see cref="IdentificationRule"/> on the value its own
        /// FDR values say (the PEP q-value only where PEP was trained for its digestion agent), strictly below the
        /// threshold. No count rule and no silent fallback.
        /// </summary>
        private static FilteredPsms FilterByName(List<SpectralMatch> candidates, CommonParameters commonParams, bool includeHighQValuePsms,
            double qValueThreshold, double pepQValueThreshold, bool filterAtPeptideLevel)
        {
            IdentificationFilterMode mode = commonParams.IdentificationFilter;
            double threshold = IdentificationRule.Threshold(mode, qValueThreshold, pepQValueThreshold);
            if (!includeHighQValuePsms)
            {
                IdentificationRule.ThrowIfPepNotTrained(candidates, filterAtPeptideLevel, commonParams);
            }

            var judgedOn = candidates
                .Where(psm => psm.GetFdrInfo(filterAtPeptideLevel) != null)
                .Select(psm => IdentificationRule.ValueFor(psm.GetFdrInfo(filterAtPeptideLevel), mode))
                .Distinct()
                .OrderBy(v => v)
                .ToList();
            FilterType filterType = judgedOn.Count switch
            {
                0 => IdentificationRule.ValueFor(null, mode),
                1 => judgedOn[0],
                _ => FilterType.QValueNotch
            };

            List<SpectralMatch> filtered = candidates
                .Where(psm => includeHighQValuePsms || IdentificationRule.Passes(psm, filterAtPeptideLevel, mode, threshold))
                .CollapseToPeptides(filterAtPeptideLevel)
                .ToList();

            return new FilteredPsms(filtered, filterType, threshold, false, filterAtPeptideLevel)
            {
                IdentificationFilter = mode,
                JudgedOn = string.Join(" or ", (judgedOn.Count == 0 ? new List<FilterType> { filterType } : judgedOn).Select(GetFilterTypeString))
            };
        }

        /// <summary>
        /// True when no match carries a computed PEP q-value at this level: each is still FdrInfo's initial 2.
        /// </summary>
        private static bool NoPepQValues(List<SpectralMatch> psms, bool filterAtPeptideLevel) =>
            psms.All(p => p.GetFdrInfo(filterAtPeptideLevel) != null && Math.Abs(p.GetFdrInfo(filterAtPeptideLevel).PEP_QValue - 2) < 0.1);

        /// <summary>
        /// What was done when PEP q-values were missing. At the peptide level PEP may still have been trained on the
        /// PSMs, so the sentence names the level, not PEP itself.
        /// </summary>
        private static string PepNotComputedSentence(bool filterAtPeptideLevel, string label) => filterAtPeptideLevel
            ? "PEP filtering was requested, but PEP q-values were not computed at the peptide level, so results were filtered by q-value."
            : $"PEP filtering was requested, but PEP was not computed for these {label}s, so results were filtered by q-value.";

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
