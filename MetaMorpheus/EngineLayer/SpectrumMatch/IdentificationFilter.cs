using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using EngineLayer.FdrAnalysis;
using Quantification;

namespace EngineLayer.SpectrumMatch
{
    /// <summary>
    /// The one identification filter of the tiered mode: the decision of which FDR value decides whether a
    /// spectral match (or peptide) counts as identified, made in one place for every consumer -- the PSM and
    /// peptide tables, protein inference and protein FDR, the PSMs handed to FlashLFQ, and the summaries.
    ///
    /// <para>Tiered mode is opt-in (<see cref="IsTieredMode(double, double)"/>): it is on when the PEP q-value
    /// threshold is below 1 and not above the q-value threshold. With the default settings (q 0.01, PEP q 1.0)
    /// it is off, and the legacy filter in <see cref="FilteredPsms"/> is used unchanged.</para>
    ///
    /// <para>In tiered mode, for each level (PSM, peptide) the first usable tier wins:
    /// <list type="number">
    /// <item>PEP q-value, if PEP was actually trained for that level (<see cref="FilterType.PepQValue"/>);</item>
    /// <item>else the q-value notch, if notch q-values were computed (<see cref="FilterType.QValueNotch"/>);</item>
    /// <item>else the q-value alone (<see cref="FilterType.QValue"/>).</item>
    /// </list>
    /// The choice depends on what was computed, never on how many matches there are, and a tier always filters.
    /// Every tier is a strict comparison: value &lt; threshold.</para>
    ///
    /// <para>The rule itself (tier order, what makes PEP or the notch usable, the strict comparison) is mzLib's
    /// <see cref="QuantifiedPsmRule"/>: <see cref="QuantifiedPsmRule.ChooseTier"/> picks the tier and
    /// <see cref="QuantifiedPsmRule.PassesConfidence"/> judges each match. This class decides what mzLib does not:
    /// whether tiered mode is on, the threshold, which level (PSM or peptide) and which matches the tier is chosen
    /// from, and the words results.txt uses to say which tier decided and why.</para>
    /// </summary>
    public static class IdentificationFilter
    {
        /// <summary>
        /// True when the tiered identification filter is in effect: the PEP q-value threshold is below 1
        /// (i.e. the user asked for PEP filtering) and not above the q-value threshold.
        /// </summary>
        public static bool IsTieredMode(double qValueThreshold, double pepQValueThreshold) =>
            pepQValueThreshold < 1.0 && pepQValueThreshold <= qValueThreshold;

        /// <inheritdoc cref="IsTieredMode(double, double)"/>
        public static bool IsTieredMode(CommonParameters commonParameters) =>
            commonParameters != null && IsTieredMode(commonParameters.QValueThreshold, commonParameters.PepQValueThreshold);

        /// <summary>
        /// The threshold every tier is held to: the smaller of the two configured thresholds. In tiered mode that
        /// is the PEP q-value threshold. A q-value threshold of 1 (what the GUI sends when its PEP box is ticked)
        /// therefore never turns a fallback tier into "filter nothing".
        /// </summary>
        public static double Threshold(double qValueThreshold, double pepQValueThreshold) =>
            Math.Min(qValueThreshold, pepQValueThreshold);

        /// <summary>
        /// Chooses the tier for one level of one set of matches, at the threshold implied by the parameters.
        /// </summary>
        public static IdentificationTier Resolve(IEnumerable<SpectralMatch> matches, bool peptideLevel, CommonParameters commonParameters) =>
            Resolve(matches, peptideLevel, Threshold(commonParameters.QValueThreshold, commonParameters.PepQValueThreshold));

        /// <summary>
        /// Chooses the tier for one level of one set of matches.
        /// </summary>
        /// <param name="matches">the matches the tier is decided for (nulls are ignored)</param>
        /// <param name="peptideLevel">true to decide on the peptide-level FDR values, false for the PSM level</param>
        /// <param name="threshold">the threshold the chosen value must be strictly below</param>
        public static IdentificationTier Resolve(IEnumerable<SpectralMatch> matches, bool peptideLevel, double threshold)
        {
            List<SpectralMatch> withFdr = (matches ?? Enumerable.Empty<SpectralMatch>())
                .Where(m => m?.GetFdrInfo(peptideLevel) != null)
                .ToList();

            // PEP q-values: targets only, so decoys alone never make PEP look trained.
            // PEP itself lives on the PSM-level FdrInfo, whichever level is being filtered.
            List<double> pepQValues = withFdr.Where(m => !m.IsDecoy).Select(m => m.GetFdrInfo(peptideLevel).PEP_QValue).ToList();
            List<double> notchQValues = withFdr.Select(m => m.GetFdrInfo(peptideLevel).QValueNotch).ToList();
            List<double> peps = withFdr.Where(m => m.PsmFdrInfo != null).Select(m => m.PsmFdrInfo.PEP).ToList();

            FilterType filterType = ToFilterType(QuantifiedPsmRule.ChooseTier(pepQValues, notchQValues, peps));
            if (filterType == FilterType.PepQValue)
            {
                return new IdentificationTier(FilterType.PepQValue, threshold, peptideLevel, null);
            }

            string pepReason = PepNotUsableReason(withFdr, pepQValues, peps, peptideLevel);
            return filterType == FilterType.QValueNotch
                ? new IdentificationTier(FilterType.QValueNotch, threshold, peptideLevel, pepReason)
                : new IdentificationTier(FilterType.QValue, threshold, peptideLevel, pepReason + "; no q-value notch was computed");
        }

        /// <summary>The MetaMorpheus filter type for mzLib's tier: each tier names the same value.</summary>
        public static FilterType ToFilterType(QuantifiedPsmTier tier) => tier switch
        {
            QuantifiedPsmTier.PepQValue => FilterType.PepQValue,
            QuantifiedPsmTier.QValueNotch => FilterType.QValueNotch,
            QuantifiedPsmTier.QValue => FilterType.QValue,
            _ => throw new ArgumentOutOfRangeException(nameof(tier), tier, null)
        };

        /// <summary>mzLib's tier for a MetaMorpheus filter type: the inverse of <see cref="ToFilterType"/>.</summary>
        public static QuantifiedPsmTier ToTier(FilterType filterType) => filterType switch
        {
            FilterType.PepQValue => QuantifiedPsmTier.PepQValue,
            FilterType.QValueNotch => QuantifiedPsmTier.QValueNotch,
            FilterType.QValue => QuantifiedPsmTier.QValue,
            _ => throw new ArgumentOutOfRangeException(nameof(filterType), filterType, null)
        };

        /// <summary>
        /// The value a tier compares against its threshold for one set of FDR values.
        /// </summary>
        public static double GetValue(FdrInfo fdrInfo, FilterType filterType) => filterType switch
        {
            FilterType.PepQValue => fdrInfo.PEP_QValue,
            FilterType.QValueNotch => fdrInfo.QValueNotch,
            _ => fdrInfo.QValue
        };

        /// <summary>
        /// The tiered rule for one set of FDR values: mzLib's <see cref="QuantifiedPsmRule.PassesConfidence"/>, the
        /// tier's value strictly below the threshold. FdrInfo starts every value at 2, so a value that was never
        /// computed never passes.
        /// </summary>
        public static bool Passes(FdrInfo fdrInfo, FilterType filterType, double threshold) =>
            fdrInfo != null
            && QuantifiedPsmRule.PassesConfidence(ToTier(filterType), fdrInfo.QValue, fdrInfo.QValueNotch, fdrInfo.PEP_QValue, threshold);

        /// <summary>
        /// Why <see cref="QuantifiedPsmRule.PepIsUsable"/> said no, in words for results.txt. Either no target
        /// carries a PEP q-value in [0, 1] (PEP was not trained), or the PEP values are all one value (a failed
        /// training run: PEPAnalysisEngine returns early, every PEP stays 0, and the FDR engine still computes a PEP
        /// q-value that is only a score ranking).
        /// </summary>
        private static string PepNotUsableReason(List<SpectralMatch> withFdr, List<double> targetPepQValues,
            List<double> peps, bool peptideLevel)
        {
            if (!QuantifiedPsmRule.PepIsUsable(targetPepQValues))
            {
                int count = peptideLevel
                    ? withFdr.Select(m => m.FullSequence).Distinct().Count()
                    : withFdr.Count;
                string unit = peptideLevel
                    ? GlobalVariables.AnalyteType.GetUniqueFormLabel().ToLower() + "s"
                    : GlobalVariables.AnalyteType.GetSpectralMatchLabel() + "s";
                return $"PEP not trained: {count} {unit}";
            }

            double onlyPep = peps.First(double.IsFinite);
            return $"PEP training failed: every PEP is {onlyPep.ToString(CultureInfo.InvariantCulture)}";
        }
    }

    /// <summary>
    /// The tier <see cref="IdentificationFilter.Resolve(IEnumerable{SpectralMatch}, bool, double)"/> chose for one level.
    /// </summary>
    public sealed class IdentificationTier
    {
        public IdentificationTier(FilterType filterType, double threshold, bool peptideLevel, string fallbackReason)
        {
            FilterType = filterType;
            Threshold = threshold;
            PeptideLevel = peptideLevel;
            FallbackReason = fallbackReason;
        }

        /// <summary>Which value decides: PEP q-value, q-value notch, or q-value.</summary>
        public FilterType FilterType { get; }

        /// <summary>The value must be strictly below this.</summary>
        public double Threshold { get; }

        /// <summary>True when the tier reads the peptide-level FDR values.</summary>
        public bool PeptideLevel { get; }

        /// <summary>Null when the PEP q-value is used; otherwise why it could not be.</summary>
        public string FallbackReason { get; }

        /// <summary>True when the match passes this tier (its value strictly below the threshold).</summary>
        public bool Passes(SpectralMatch match) =>
            match != null && IdentificationFilter.Passes(match.GetFdrInfo(PeptideLevel), FilterType, Threshold);

        /// <summary>The value this tier reads for a match (for example, the q-value handed to FlashLFQ).</summary>
        public double GetValue(SpectralMatch match) => IdentificationFilter.GetValue(match.GetFdrInfo(PeptideLevel), FilterType);

        /// <summary>E.g. "pep q-value &lt; 0.01".</summary>
        public string Describe() => $"{FilteredPsms.GetFilterTypeString(FilterType)} < {Threshold.ToString(CultureInfo.InvariantCulture)}";

        /// <summary>E.g. "q-value notch &lt; 0.01 (PEP not trained: 640 PSMs)".</summary>
        public string DescribeWithReason() => FallbackReason == null ? Describe() : $"{Describe()} ({FallbackReason})";
    }
}
