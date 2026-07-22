using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace EngineLayer.Deconvolution.FeatureFileMapping;

/// <summary>
/// One raw spectra file and all of its condition-specific feature-file associations.
/// </summary>
public class FeatureSpectraEntry
{
    public string MassSpecFilePath { get; set; } = string.Empty;

    public string MassSpecFileName { get; set; } = string.Empty;

    public string MassSpecFileNameWithoutExtension { get; set; } = string.Empty;

    public long? MassSpecFileSizeBytes { get; set; }

    /// <summary>
    /// All feature-file variants available for this raw file, keyed by condition.
    /// </summary>
    public List<FeatureFileConditionEntry> FeatureFiles { get; set; } = new();

    public static FeatureSpectraEntry Create(string rawFilePath)
    {
        string normalizedPath = NormalizePath(rawFilePath);
        return new FeatureSpectraEntry
        {
            MassSpecFilePath = normalizedPath,
            MassSpecFileName = Path.GetFileName(normalizedPath),
            MassSpecFileNameWithoutExtension = Path.GetFileNameWithoutExtension(normalizedPath),
            MassSpecFileSizeBytes = File.Exists(normalizedPath) ? new FileInfo(normalizedPath).Length : null
        };
    }

    public bool MatchesRawFile(string rawFilePath)
    {
        if (string.IsNullOrEmpty(rawFilePath))
            return string.IsNullOrEmpty(MassSpecFilePath);
        return string.Equals(MassSpecFilePath, NormalizePath(rawFilePath), StringComparison.OrdinalIgnoreCase);
    }

    private static string NormalizePath(string path)
    {
        if (string.IsNullOrEmpty(path))
            return path ?? string.Empty;
        return Path.GetFullPath(path);
    }

    public void AddOrReplaceConditionFile(FeatureFileConditionEntry conditionEntry)
    {
        ArgumentNullException.ThrowIfNull(conditionEntry);

        int existingIndex = FeatureFiles.FindIndex(p => p.MatchesCondition(conditionEntry.ConditionKey));
        if (existingIndex >= 0)
        {
            FeatureFiles[existingIndex] = conditionEntry;
        }
        else
        {
            FeatureFiles.Add(conditionEntry);
        }
    }

    public bool TryGetConditionFile(string conditionKey, out FeatureFileConditionEntry entry)
    {
        entry = FeatureFiles.FirstOrDefault(p => p.MatchesCondition(conditionKey));
        return entry is not null;
    }
}
