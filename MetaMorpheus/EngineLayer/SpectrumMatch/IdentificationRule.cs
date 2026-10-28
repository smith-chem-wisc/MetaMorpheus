using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using EngineLayer.FdrAnalysis;
using Quantification;

namespace EngineLayer.SpectrumMatch
{
    /// <summary>
    /// The identification filter a search asked for by name (<see cref="CommonParameters.IdentificationFilter"/>
    /// other than <see cref="IdentificationFilterMode.Thresholds"/>), in one place for every consumer: the PSM and
    /// peptide tables, protein inference and protein FDR, quantification, the pruned database, the spectral library,
    /// the mzIdentML and pepXML files, and the summaries.
    ///
    /// <para>The FDR value each match is judged on is decided by the match's own <see cref="FdrInfo.PepTrained"/>,
    /// which <see cref="FdrAnalysisEngine"/> sets for each digestion agent at each level. So the decision is made
    /// once, by the engine that trained (or did not train) PEP, and every later filter over any subset of the
    /// matches reaches the same answer for each match. A digestion agent PEP was not trained for is never judged on
    /// a PEP q-value.</para>
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
        /// PEP q-value, and for PEP-else-notch the smaller of the two, for both tiers. A q-value threshold of 1 (what
        /// the GUI sends when its PEP box is ticked) therefore never turns the fallback into "filter nothing".
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

        /// <summary>The FDR value a match is judged on, under <paramref name="mode"/>, given its FDR values at one level.</summary>
        public static FilterType ValueFor(FdrInfo fdrInfo, IdentificationFilterMode mode) => mode switch
        {
            IdentificationFilterMode.QValueNotch => FilterType.QValueNotch,
            IdentificationFilterMode.PepQValue => FilterType.PepQValue,
            IdentificationFilterMode.PepQValueElseQValueNotch => fdrInfo is { PepTrained: true } ? FilterType.PepQValue : FilterType.QValueNotch,
            _ => throw new ArgumentOutOfRangeException(nameof(mode), mode, "Not an explicit identification filter.")
        };

        /// <summary>
        /// True when the match passes the filter at this level: its value strictly below <paramref name="threshold"/>.
        /// A PEP q-value that did not come from a trained model never passes.
        /// </summary>
        public static bool Passes(SpectralMatch match, bool peptideLevel, IdentificationFilterMode mode, double threshold)
        {
            FdrInfo fdrInfo = match?.GetFdrInfo(peptideLevel);
            if (fdrInfo == null)
            {
                return false;
            }
            FilterType value = ValueFor(fdrInfo, mode);
            if (value == FilterType.PepQValue && !fdrInfo.PepTrained)
            {
                return false;
            }
            return QuantifiedPsmRule.PassesConfidence(ToTier(value), fdrInfo.QValue, fdrInfo.QValueNotch, fdrInfo.PEP_QValue, threshold);
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
        /// For each digestion agent among <paramref name="matches"/>, whether PEP was trained at this level, and how
        /// many matches (PSM level) or distinct full sequences (peptide level) it has. In name order.
        /// </summary>
        public static List<(string Agent, bool PepTrained, int Count)> PepTrainingByAgent(IEnumerable<SpectralMatch> matches, bool peptideLevel) =>
            matches
                .Where(m => m?.GetFdrInfo(peptideLevel) != null)
                .GroupBy(m => m.DigestionParams?.DigestionAgent?.Name ?? "unknown")
                .OrderBy(g => g.Key, StringComparer.Ordinal)
                .Select(g => (g.Key,
                    g.Any(m => m.GetFdrInfo(peptideLevel).PepTrained),
                    peptideLevel ? g.Select(m => m.FullSequence).Distinct().Count() : g.Count()))
                .ToList();

        /// <summary>
        /// Under <see cref="IdentificationFilterMode.PepQValue"/>, the search cannot be filtered as asked when PEP was
        /// not trained for a digestion agent: this throws, naming each such agent. Does nothing under any other mode.
        /// </summary>
        public static void ThrowIfPepNotTrained(IEnumerable<SpectralMatch> matches, bool peptideLevel, CommonParameters commonParameters)
        {
            if (commonParameters?.IdentificationFilter != IdentificationFilterMode.PepQValue)
            {
                return;
            }
            var untrained = PepTrainingByAgent(matches, peptideLevel).Where(a => !a.PepTrained).ToList();
            if (untrained.Count > 0)
            {
                throw new MetaMorpheusException(
                    $"The identification filter is the PEP q-value, but PEP was not trained at the {LevelName(peptideLevel)} level for "
                    + string.Join(", ", untrained.Select(a => $"{a.Agent} ({a.Count} {UnitName(peptideLevel)})"))
                    + ". PEP is trained only for a digestion agent with more than 1000 matches, and training can fail on small or one-sided data. "
                    + "Choose PepQValueElseQValueNotch or QValueNotch as the identification filter to search these data.");
            }
        }

        /// <summary>
        /// What the filter did at this level, for results.txt and the manuscript prose: e.g. "pep q-value &lt; 0.01",
        /// or, when digestion agents differ, "pep q-value &lt; 0.01 for trypsin; q-value notch &lt; 0.01 for Arg-C (PEP
        /// not trained: 640 PSMs)".
        /// </summary>
        public static string Describe(IEnumerable<SpectralMatch> matches, bool peptideLevel, CommonParameters commonParameters)
        {
            IdentificationFilterMode mode = commonParameters.IdentificationFilter;
            string threshold = " < " + Threshold(commonParameters).ToString(CultureInfo.InvariantCulture);
            if (mode != IdentificationFilterMode.PepQValueElseQValueNotch)
            {
                return FilteredPsms.GetFilterTypeString(ValueFor(null, mode)) + threshold;
            }

            var agents = PepTrainingByAgent(matches, peptideLevel);
            if (agents.All(a => a.PepTrained))
            {
                return FilteredPsms.GetFilterTypeString(FilterType.PepQValue) + threshold;
            }
            string Fallback((string Agent, bool PepTrained, int Count) a) =>
                $"{FilteredPsms.GetFilterTypeString(FilterType.QValueNotch)}{threshold} (PEP not trained: {a.Count} {UnitName(peptideLevel)})";
            if (agents.Count == 1)
            {
                return Fallback(agents[0]);
            }
            return string.Join("; ", agents.Select(a => (a.PepTrained ? FilteredPsms.GetFilterTypeString(FilterType.PepQValue) + threshold : Fallback(a)) + " for " + a.Agent));
        }

        /// <summary>True when, under PEP-else-notch, at least one digestion agent fell back to the notch at this level.</summary>
        public static bool FellBack(IEnumerable<SpectralMatch> matches, bool peptideLevel, CommonParameters commonParameters) =>
            commonParameters.IdentificationFilter == IdentificationFilterMode.PepQValueElseQValueNotch
            && PepTrainingByAgent(matches, peptideLevel).Any(a => !a.PepTrained);

        private static string LevelName(bool peptideLevel) =>
            peptideLevel ? GlobalVariables.AnalyteType.GetUniqueFormLabel().ToLower() : GlobalVariables.AnalyteType.GetSpectralMatchLabel();

        private static string UnitName(bool peptideLevel) => LevelName(peptideLevel) + "s";
    }
}
