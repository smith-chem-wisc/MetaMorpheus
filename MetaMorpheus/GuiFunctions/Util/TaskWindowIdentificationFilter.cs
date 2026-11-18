using System.Globalization;
using EngineLayer;
using EngineLayer.SpectrumMatch;

namespace GuiFunctions.Util;

/// <summary>
/// What the Search task window shows and saves for <see cref="CommonParameters.IdentificationFilter"/>: its PEP box is the
/// filter (ticked, the PEP q-value; unticked, the q-value notch), with one threshold box for each.
/// </summary>
/// <remarks>
/// A task saved before the filter was named, with the old PEP box ticked, names no filter and holds QValueThreshold 1.0 and
/// PepQValueThreshold 0.01. It loads as the PEP q-value filter (<see cref="CommonParameters.IdentificationFilterFromLegacyThresholds"/>),
/// so the window shows the box ticked, and saving it writes IdentificationFilter = "PepQValue".
/// </remarks>
public static class TaskWindowIdentificationFilter
{
    /// <summary>Whether the window shows its PEP box ticked for this loaded task.</summary>
    public static bool PepBoxTicked(CommonParameters loaded) => loaded.IdentificationFilter == FilterType.PepQValue;

    /// <summary>The q-value notch threshold the window shows: the task's, or 0.01 when the task holds 1 or more.</summary>
    public static string QValueThresholdText(CommonParameters loaded) =>
        loaded.QValueThreshold < 1 ? loaded.QValueThreshold.ToString(CultureInfo.InvariantCulture) : "0.01";

    /// <summary>The PEP q-value threshold the window shows: the task's, or 0.01 when the task holds 1 or more.</summary>
    public static string PepQValueThresholdText(CommonParameters loaded) =>
        loaded.PepQValueThreshold < 1 ? loaded.PepQValueThreshold.ToString(CultureInfo.InvariantCulture) : "0.01";

    /// <summary>
    /// The filter and the two thresholds the window saves: the box's filter at its own threshold, and 1.0 for the other,
    /// as the window has always saved them.
    /// </summary>
    public static (FilterType IdentificationFilter, double QValueThreshold, double PepQValueThreshold) ForSave(
        bool pepBoxTicked, string qValueThresholdText, string pepQValueThresholdText) =>
        pepBoxTicked
            ? (FilterType.PepQValue, 1.0, double.Parse(pepQValueThresholdText, CultureInfo.InvariantCulture))
            : (FilterType.QValueNotch, double.Parse(qValueThresholdText, CultureInfo.InvariantCulture), 1.0);
}
