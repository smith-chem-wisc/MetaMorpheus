#nullable enable
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.SpectrumMatch;

/// <summary>
/// Decides which of one match's ambiguous hypotheses PEP has scored far enough below the best to remove.
/// <remarks>
/// Run by <see cref="DisambiguationEngine"/> after FDR, on hypotheses whose <see cref="SpectralMatchHypothesis.PEP"/>
/// the PEP engine has set. Never returns the best hypothesis, so the match's PEP is unchanged.
/// </remarks>
/// </summary>
public interface IPepDisambiguationRule
{
    /// <param name="hypotheses">Every hypothesis of one match, each with a PEP.</param>
    /// <returns>The hypotheses to remove.</returns>
    IEnumerable<SpectralMatchHypothesis> HypothesesToRemove(IReadOnlyList<SpectralMatchHypothesis> hypotheses);
}

/// <summary>
/// Removes every hypothesis whose PEP is more than <see cref="MaxPepGap"/> above the best one's.
/// </summary>
public sealed class AbsolutePepGapRule(double maxPepGap) : IPepDisambiguationRule
{
    /// <summary>
    /// The rule the PEP engine applied itself before it stopped removing hypotheses: a gap of 0.05. The value has no
    /// recorded rationale. Glyco, crosslink and nonspecific searches keep it until a rule is chosen on data.
    /// </summary>
    public static AbsolutePepGapRule PepEngineRule { get; } = new(0.05);

    public double MaxPepGap { get; } = maxPepGap;

    public IEnumerable<SpectralMatchHypothesis> HypothesesToRemove(IReadOnlyList<SpectralMatchHypothesis> hypotheses)
    {
        double bestPep = hypotheses.Min(h => h.PEP!.Value);
        return hypotheses.Where(h => h.PEP!.Value - bestPep > MaxPepGap);
    }
}
