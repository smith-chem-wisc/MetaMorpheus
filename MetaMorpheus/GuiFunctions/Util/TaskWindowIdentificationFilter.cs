using System.Collections.Generic;
using System.Globalization;
using EngineLayer;
using EngineLayer.SpectrumMatch;

namespace GuiFunctions.Util;

/// <summary>
/// What a task window shows and saves for <see cref="CommonParameters.IdentificationFilter"/> and
    /// <see cref="CommonParameters.FilterThreshold"/>: a choice of the filters (<see cref="Choices"/>) and one threshold
/// box. The Search window loads them with <see cref="Load(CommonParameters, bool)"/> and saves what it shows with
/// <see cref="ForSave"/>.
/// </summary>
/// <remarks>
/// <para>A task saved before FilterThreshold existed is read by its old thresholds (an old GUI PEP task, QValueThreshold 1.0
/// and PepQValueThreshold 0.01, as the PEP q-value at 0.01), so the window shows what the task searches with, and saving it
/// writes IdentificationFilter and FilterThreshold.</para>
/// </remarks>
public static class TaskWindowIdentificationFilter
{
    /// <summary>One identification filter a task window offers: its filter, its name in the window, and what it means.</summary>
    public sealed record Choice(FilterType Filter, string Label, string Description);

    /// <summary>The filters, in the order the windows list them.</summary>
    public static IReadOnlyList<Choice> Choices { get; } = new[]
    {
        new Choice(FilterType.PepQValue, "PEP q-value",
            "The PEP q-value: matches ranked by their posterior error probability (PEP), decoys counted down that ranking. " +
            "Where PEP could not be trained for a match (too few matches, or a protease or cleavage category without enough), " +
            "the match is judged by its q-value notch at the same threshold instead, and the search says so. The default bottom-up."),
        new Choice(FilterType.QValueNotch, "q-value notch",
            "The q-value counted within each mass-difference notch (each precursor mass offset the search accepts). " +
            "The default for top-down and RNA."),
        new Choice(FilterType.QValue, "q-value (pooled)",
            "The q-value counted over every mass-difference notch together, and nothing else."),
        new Choice(FilterType.QValueAndQValueNotch, "q-value and q-value notch",
            "The pooled q-value and q-value notch, both held to the same threshold."),
    };

    /// <summary>
    /// The filter and threshold a window shows for the loaded task. A task still at the bottom-up default (the PEP q-value at
    /// 0.01) shows the default for its search type instead, as the Search window does for the PEP training cutoff: the
    /// q-value notch for a top-down protease (<paramref name="isTopDown"/>) and in RNA mode (<paramref name="isRna"/>).
    /// </summary>
    public static (FilterType Filter, string ThresholdText) Load(CommonParameters loaded, bool isRna, bool isTopDown)
    {
        FilterType filter = loaded.IdentificationFilter;
        double threshold = loaded.FilterThreshold;
        if (filter == CommonParameters.DefaultIdentificationFilter(false, false) && threshold == CommonParameters.DefaultFilterThreshold)
        {
            filter = CommonParameters.DefaultIdentificationFilter(isRna, isTopDown);
        }
        return (filter, threshold.ToString(CultureInfo.InvariantCulture));
    }

    /// <summary>
    /// <see cref="Load(CommonParameters, bool, bool)"/> for the loaded task's own protease: top-down when it is "top-down".
    /// What a task window calls.
    /// </summary>
    public static (FilterType Filter, string ThresholdText) Load(CommonParameters loaded, bool isRna) =>
        Load(loaded, isRna, loaded.DigestionParams?.DigestionAgent?.Name == "top-down");

    /// <summary>
    /// The filter and threshold a window sets when the user takes the top-down recommended settings: the top-down default
    /// (<see cref="CommonParameters.DefaultIdentificationFilter"/>) at <see cref="CommonParameters.DefaultFilterThreshold"/>.
    /// </summary>
    public static (FilterType Filter, string ThresholdText) TopDownRecommended(bool isRna) =>
        (CommonParameters.DefaultIdentificationFilter(isRna, isTopDown: true),
            CommonParameters.DefaultFilterThreshold.ToString(CultureInfo.InvariantCulture));

    /// <summary>The line the top-down recommendation prompt adds for the filter, naming the filter and threshold it sets.</summary>
    public static string TopDownRecommendedSentence(bool isRna)
    {
        var (filter, thresholdText) = TopDownRecommended(isRna);
        return $"Use '{LabelOf(filter)}' at '{thresholdText}' for 'Identification Filter'";
    }

    /// <summary>The filter and threshold the window saves: what its controls show.</summary>
    public static (FilterType IdentificationFilter, double FilterThreshold) ForSave(FilterType filter, string thresholdText) =>
        (filter, double.Parse(thresholdText, NumberStyles.Float, CultureInfo.InvariantCulture));

    /// <summary>
    /// Why <paramref name="thresholdText"/> is not a threshold a filter can be held to (a number above 0 and at most 1), or
    /// null when it is one.
    /// </summary>
    public static string ThresholdError(string thresholdText) =>
        double.TryParse(thresholdText, NumberStyles.Float, CultureInfo.InvariantCulture, out double threshold) && threshold > 0 && threshold <= 1
            ? null
            : $"The identification filter threshold must be a number above 0 and at most 1. \n You entered \"{thresholdText}\"";

    /// <summary>
    /// The warning for a threshold of 1, which keeps every match (allowed: the search says so too), or null for any other
    /// threshold.
    /// </summary>
    public static string ThresholdWarning(string thresholdText) =>
        double.TryParse(thresholdText, NumberStyles.Float, CultureInfo.InvariantCulture, out double threshold) && threshold >= 1
            ? "An identification filter threshold of 1 keeps every match: nothing is filtered out of the results, protein inference, quantification or the counts."
            : null;

    private static string LabelOf(FilterType filter)
    {
        foreach (Choice choice in Choices)
        {
            if (choice.Filter == filter)
            {
                return choice.Label;
            }
        }
        return filter.ToString();
    }
}
