using System;

namespace EngineLayer.Deconvolution.FeatureFileMapping;

/// <summary>
/// Declares one named condition that can be applied across many spectra files. e.g. FlashDeconv vs Dinosaur, or Dinosaur with tight RT vs Dinosaur with loose RT.
/// </summary>
public class FeatureFileMapCondition
{
    /// <summary>
    /// Stable machine-readable key, e.g. <c>flashdeconv-default</c> or <c>dinosaur-tight-rt</c>.
    /// </summary>
    public string ConditionKey { get; set; } = string.Empty;

    /// <summary>
    /// Human-readable display name shown in UI/CLI selection surfaces.
    /// </summary>
    public string DisplayName { get; set; } = string.Empty;

    public string Description { get; set; } = string.Empty;

    public string FeatureGenerationMethod { get; set; } = string.Empty;

    public string FeatureGenerationVersion { get; set; } = string.Empty;

    public bool MatchesKey(string conditionKey)
        => string.Equals(ConditionKey, conditionKey, StringComparison.OrdinalIgnoreCase);
}
