using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using EngineLayer.FdrAnalysis;
using Omics.Digestion;
using Quantification;

namespace EngineLayer.SpectrumMatch
{
    /// <summary>
    /// The identification filter a search asked for by name (<see cref="CommonParameters.IdentificationFilter"/>
    /// other than <see cref="IdentificationFilterMode.Thresholds"/>), in one place for every consumer: the PSM and
    /// peptide tables, protein inference and protein FDR, quantification, the pruned database, the spectral library,
    /// the mzIdentML and pepXML files, and the summaries.
    ///
    /// <para>The FDR value each match is judged on is decided by the match's own <see cref="FdrInfo.PepTrained"/>.
    /// <see cref="FdrAnalysisEngine"/> sets it for each digestion agent at each level, on the matches it is given; a
    /// semi- or non-specific search gives it each cleavage-specificity category separately. So the decision is made
    /// once, by the engine that trained (or did not train) PEP, and every later filter over any subset of the
    /// matches reaches the same answer for each match. A match whose PEP was not trained is never judged on a PEP
    /// q-value: both PEP filters fall back for it (<see cref="Fallback"/>), and the search reports each fallback
    /// (<see cref="PepNotTrainedSentence"/>).</para>
    ///
    /// <para>The comparison is mzLib's <see cref="QuantifiedPsmRule.PassesConfidence"/>, the rule FlashLFQ also
    /// applies: the value strictly below the threshold.</para>
    /// </summary>
    public static class IdentificationRule
    {
        /// <summary>True when the search names its identification filter instead of inferring it from the thresholds.</summary>
        public static bool IsExplicit(CommonParameters commonParameters) =>
            commonParameters != null && commonParameters.IdentificationFilter != IdentificationFilterMode.Thresholds;

        /// <summary>
        /// The threshold the filter is held to: the q-value threshold for the notch, the PEP q-value threshold for the
        /// PEP q-value (and its fallback), and for PEP-else-notch the smaller of the two, for both tiers.
        /// </summary>
        public static double Threshold(CommonParameters commonParameters) =>
            Threshold(commonParameters.IdentificationFilter, commonParameters.QValueThreshold, commonParameters.PepQValueThreshold);

        /// <inheritdoc cref="Threshold(CommonParameters)"/>
        public static double Threshold(IdentificationFilterMode mode, double qValueThreshold, double pepQValueThreshold) => mode switch
        {
            IdentificationFilterMode.QValueNotch => qValueThreshold,
            IdentificationFilterMode.PepQValue => pepQValueThreshold,
            IdentificationFilterMode.PepQValueElseQValueNotch => Math.Min(qValueThreshold, pepQValueThreshold),
            _ => throw new ArgumentOutOfRangeException(nameof(mode), mode,
                "The identification filter is inferred from the thresholds; there is no single threshold.")
        };

        /// <summary>
        /// Why a named filter cannot run with these thresholds, or null when it can (and when no filter is named). Its
        /// threshold must lie strictly between 0 and 1: compared strictly, a threshold of 1 keeps every match, and 1 is
        /// what the GUI sends for the threshold box it leaves unticked; 0 keeps none.
        /// </summary>
        public static string ThresholdProblem(CommonParameters commonParameters)
        {
            if (!IsExplicit(commonParameters))
            {
                return null;
            }
            IdentificationFilterMode mode = commonParameters.IdentificationFilter;
            double threshold = Threshold(commonParameters);
            if (threshold > 0 && threshold < 1)
            {
                return null;
            }
            string heldTo = mode switch
            {
                IdentificationFilterMode.QValueNotch => "the q-value threshold",
                IdentificationFilterMode.PepQValue => "the PEP q-value threshold",
                _ => "the smaller of the q-value and PEP q-value thresholds"
            };
            string set = mode switch
            {
                IdentificationFilterMode.QValueNotch => "Set the q-value threshold between 0 and 1 (in the Search task window, leave the PEP q-value box unticked)",
                IdentificationFilterMode.PepQValue => "Set the PEP q-value threshold between 0 and 1 (in the Search task window, tick the PEP q-value box)",
                _ => "Set either threshold between 0 and 1"
            };
            return $"The identification filter {mode} holds every match to {heldTo}, which is {threshold.ToString(CultureInfo.InvariantCulture)}, "
                + $"so it would keep {(threshold >= 1 ? "every match" : "no match")}. {set}, or choose another identification filter.";
        }

        /// <summary>
        /// Refuses a named filter whose threshold is not strictly between 0 and 1 (<see cref="ThresholdProblem"/>).
        /// Called when a task starts, before anything is searched, so nothing is lost.
        /// </summary>
        public static void ThrowIfThresholdUnusable(CommonParameters commonParameters)
        {
            if (ThresholdProblem(commonParameters) is { } problem)
            {
                throw new MetaMorpheusException(problem);
            }
        }

        /// <summary>The FDR value a match is judged on, under <paramref name="mode"/>, given its FDR values at one level.</summary>
        public static FilterType ValueFor(FdrInfo fdrInfo, IdentificationFilterMode mode) => mode switch
        {
            IdentificationFilterMode.QValueNotch => FilterType.QValueNotch,
            IdentificationFilterMode.PepQValue => fdrInfo is { PepTrained: true } ? FilterType.PepQValue : FilterType.QValue,
            IdentificationFilterMode.PepQValueElseQValueNotch => fdrInfo is { PepTrained: true } ? FilterType.PepQValue : FilterType.QValueNotch,
            _ => throw new ArgumentOutOfRangeException(nameof(mode), mode, "Not an explicit identification filter.")
        };

        /// <summary>
        /// The value a match whose PEP was not trained is judged on: the q-value under
        /// <see cref="IdentificationFilterMode.PepQValue"/>, the q-value notch otherwise.
        /// </summary>
        public static FilterType Fallback(IdentificationFilterMode mode) => ValueFor(null, mode);

        /// <summary>
        /// True when the match passes the filter at this level: its value strictly below <paramref name="threshold"/>.
        /// A q-value is read as MetaMorpheus's q-value filter reads it, the q-value and the q-value notch together.
        /// </summary>
        public static bool Passes(SpectralMatch match, bool peptideLevel, IdentificationFilterMode mode, double threshold)
        {
            FdrInfo fdrInfo = match?.GetFdrInfo(peptideLevel);
            if (fdrInfo == null)
            {
                return false;
            }
            bool Below(FilterType value) =>
                QuantifiedPsmRule.PassesConfidence(ToTier(value), fdrInfo.QValue, fdrInfo.QValueNotch, fdrInfo.PEP_QValue, threshold);
            FilterType judgedOn = ValueFor(fdrInfo, mode);
            return judgedOn == FilterType.QValue
                ? Below(FilterType.QValue) && Below(FilterType.QValueNotch)
                : Below(judgedOn);
        }

        /// <summary>The value the filter judges the match on at this level, e.g. the q-value handed to FlashLFQ.</summary>
        public static double GetValue(SpectralMatch match, bool peptideLevel, IdentificationFilterMode mode)
        {
            FdrInfo fdrInfo = match.GetFdrInfo(peptideLevel);
            return ValueFor(fdrInfo, mode) switch
            {
                FilterType.PepQValue => fdrInfo.PEP_QValue,
                FilterType.QValueNotch => fdrInfo.QValueNotch,
                _ => fdrInfo.QValue
            };
        }

        /// <summary>mzLib's tier for a MetaMorpheus filter type: each names the same value.</summary>
        public static QuantifiedPsmTier ToTier(FilterType filterType) => filterType switch
        {
            FilterType.PepQValue => QuantifiedPsmTier.PepQValue,
            FilterType.QValueNotch => QuantifiedPsmTier.QValueNotch,
            FilterType.QValue => QuantifiedPsmTier.QValue,
            _ => throw new ArgumentOutOfRangeException(nameof(filterType), filterType, null)
        };

        /// <summary>
        /// Matches PEP was trained, or not trained, for together. <see cref="Category"/> is the cleavage-specificity
        /// category, named only when the digestion agent's categories differ, which happens when FDR runs per category
        /// (semi- and non-specific searches); null otherwise. <see cref="Count"/> is matches at the PSM level and
        /// distinct full sequences at the peptide level.
        /// </summary>
        public readonly record struct PepTrainingGroup(string Agent, string Category, bool PepTrained, int Count)
        {
            /// <summary>The digestion agent, and the category when there is one, e.g. "trypsin semi-specific".</summary>
            public string Name => Category == null ? Agent : Agent + " " + Category;
        }

        /// <summary>
        /// Whether PEP was trained at this level, judged on each match's own record: one group per digestion agent
        /// when all its matches agree, and one per cleavage-specificity category when they do not. Agents in name
        /// order, categories fully, semi-, then non-specific.
        /// </summary>
        public static List<PepTrainingGroup> PepTrainingGroups(IEnumerable<SpectralMatch> matches, bool peptideLevel)
        {
            int Count(IEnumerable<SpectralMatch> group) =>
                peptideLevel ? group.Select(m => m.FullSequence).Distinct().Count() : group.Count();
            bool Trained(SpectralMatch m) => m.GetFdrInfo(peptideLevel).PepTrained;

            var groups = new List<PepTrainingGroup>();
            var byAgent = matches
                .Where(m => m?.GetFdrInfo(peptideLevel) != null)
                .GroupBy(m => m.DigestionParams.DigestionAgent.Name) // as FdrAnalysisEngine groups them
                .OrderBy(g => g.Key, StringComparer.Ordinal);
            foreach (var agent in byAgent)
            {
                var flags = agent.Select(Trained).Distinct().ToList();
                if (flags.Count == 1)
                {
                    groups.Add(new PepTrainingGroup(agent.Key, null, flags[0], Count(agent)));
                    continue;
                }
                foreach (var category in agent.GroupBy(m => (Category: CategoryOf(m), Trained: Trained(m)))
                             .OrderBy(g => g.Key.Category.Rank).ThenBy(g => g.Key.Category.Label, StringComparer.Ordinal).ThenByDescending(g => g.Key.Trained))
                {
                    groups.Add(new PepTrainingGroup(agent.Key, category.Key.Category.Label, category.Key.Trained, Count(category)));
                }
            }
            return groups;
        }

        /// <summary>
        /// What the filter did at this level, for results.txt and the manuscript prose: e.g. "pep q-value &lt; 0.01",
        /// or, when the groups differ, "pep q-value &lt; 0.01 for trypsin; q-value &lt; 0.01 (PEP not trained: 640 PSMs)
        /// for Arg-C".
        /// </summary>
        public static string Describe(IEnumerable<SpectralMatch> matches, bool peptideLevel, CommonParameters commonParameters)
        {
            IdentificationFilterMode mode = commonParameters.IdentificationFilter;
            string threshold = " < " + Threshold(commonParameters).ToString(CultureInfo.InvariantCulture);
            if (mode == IdentificationFilterMode.QValueNotch)
            {
                return FilteredPsms.GetFilterTypeString(FilterType.QValueNotch) + threshold;
            }

            var groups = PepTrainingGroups(matches, peptideLevel);
            if (groups.All(g => g.PepTrained))
            {
                return FilteredPsms.GetFilterTypeString(FilterType.PepQValue) + threshold;
            }
            string Untrained(PepTrainingGroup g) =>
                $"{FilteredPsms.GetFilterTypeString(Fallback(mode))}{threshold} (PEP not trained: {g.Count} {UnitName(peptideLevel)})";
            if (groups.Count == 1)
            {
                return Untrained(groups[0]);
            }
            return string.Join("; ", groups.Select(g => (g.PepTrained ? FilteredPsms.GetFilterTypeString(FilterType.PepQValue) + threshold : Untrained(g)) + " for " + g.Name));
        }

        /// <summary>True when, under a PEP filter, the matches of at least one group fell back at this level.</summary>
        public static bool FellBack(IEnumerable<SpectralMatch> matches, bool peptideLevel, CommonParameters commonParameters) =>
            IsPepFilter(commonParameters) && PepTrainingGroups(matches, peptideLevel).Any(g => !g.PepTrained);

        /// <summary>
        /// Under a PEP filter, the sentence saying which matches PEP was not trained for and what they were filtered by
        /// instead, in the words a PEP filter that falls back uses everywhere (<see cref="FilteredPsms.PepFilteringNotDone"/>),
        /// naming each digestion agent, and its category where there is one. Null when nothing fell back.
        /// </summary>
        public static string PepNotTrainedSentence(IEnumerable<SpectralMatch> matches, bool peptideLevel, CommonParameters commonParameters)
        {
            if (!IsPepFilter(commonParameters))
            {
                return null;
            }
            var untrained = PepTrainingGroups(matches, peptideLevel).Where(g => !g.PepTrained).Select(g => g.Name).ToList();
            if (untrained.Count == 0)
            {
                return null;
            }
            string names = untrained.Count == 1
                ? untrained[0]
                : string.Join(", ", untrained.Take(untrained.Count - 1)) + " and " + untrained[^1];
            string fallback = FilteredPsms.GetFilterTypeString(Fallback(commonParameters.IdentificationFilter));
            return peptideLevel
                ? $"PEP filtering was requested, but PEP q-values were not computed at the peptide level for {names}, so those {UnitName(true)} were filtered by {fallback}."
                : $"PEP filtering was requested, but PEP was not computed for the {names} {UnitName(false)}, so they were filtered by {fallback}.";
        }

        private static bool IsPepFilter(CommonParameters commonParameters) =>
            commonParameters?.IdentificationFilter is IdentificationFilterMode.PepQValue or IdentificationFilterMode.PepQValueElseQValueNotch;

        /// <summary>The match's cleavage-specificity category, from its best candidate (a match always has one).</summary>
        private static (int Rank, string Label) CategoryOf(SpectralMatch match) =>
            match.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer.CleavageSpecificityForFdrCategory switch
            {
                CleavageSpecificity.Full => (0, "fully specific"),
                CleavageSpecificity.Semi => (1, "semi-specific"),
                CleavageSpecificity.None => (2, "non-specific"),
                var other => (3, other.ToString())
            };

        private static string UnitName(bool peptideLevel) =>
            (peptideLevel ? GlobalVariables.AnalyteType.GetUniqueFormLabel().ToLower() : GlobalVariables.AnalyteType.GetSpectralMatchLabel()) + "s";
    }
}
