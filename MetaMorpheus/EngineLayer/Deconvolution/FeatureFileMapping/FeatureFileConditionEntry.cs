using System;
using System.IO;

namespace EngineLayer.Deconvolution.FeatureFileMapping;

/// <summary>
/// A condition-specific feature-file association for a single raw spectra file.
/// </summary>
public class FeatureFileConditionEntry
{
    public string ConditionKey { get; set; } = string.Empty;

    public string FeatureFilePath { get; set; } = string.Empty;

    public string FeatureFileName { get; set; } = string.Empty;

    public long? FeatureFileSizeBytes { get; set; }

    public static FeatureFileConditionEntry Create(string conditionKey, string featureFilePath)
    {
        string normalizedPath = FeatureFileMapStore.NormalizePath(featureFilePath);
        return new FeatureFileConditionEntry
        {
            ConditionKey = conditionKey,
            FeatureFilePath = normalizedPath,
            FeatureFileName = Path.GetFileName(normalizedPath),
            FeatureFileSizeBytes = File.Exists(normalizedPath) ? new FileInfo(normalizedPath).Length : null
        };
    }

    public bool MatchesCondition(string conditionKey)
        => string.Equals(ConditionKey, conditionKey, StringComparison.OrdinalIgnoreCase);
}
