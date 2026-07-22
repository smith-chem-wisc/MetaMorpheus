using Nett;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Deconvolution.FeatureFileMapping;

/// <summary>
/// A persisted raw-to-feature association document. One file can contain many raw
/// files and many feature-file conditions, so the same spectra can be searched
/// against multiple feature sets without creating raw-adjacent TOMLs.
/// </summary>
public class FeatureFileMap : IEnumerable<FeatureSpectraEntry>
{
    [TomlIgnore] public string FilePath { get; set; } = string.Empty;

    public DateTime LastModifiedUtc { get; set; } = DateTime.UtcNow;

    /// <summary>
    /// The named conditions available in this map, e.g. FlashDeconv vs Dinosaur.
    /// </summary>
    public List<FeatureFileMapCondition> Conditions { get; set; } = new();

    /// <summary>
    /// One entry per spectra file, each with zero or more condition-specific feature files.
    /// </summary>
    public List<FeatureSpectraEntry> SpectraFiles { get; set; } = new();

    [TomlIgnore]
    public int Count => SpectraFiles.Count;

    public IEnumerator<FeatureSpectraEntry> GetEnumerator() => SpectraFiles.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    #region Dictionary-like Methods

    public void AddOrReplaceCondition(FeatureFileMapCondition condition)
    {
        ArgumentNullException.ThrowIfNull(condition);

        int existingIndex = Conditions.FindIndex(p => p.MatchesKey(condition.ConditionKey));
        if (existingIndex >= 0)
        {
            Conditions[existingIndex] = condition;
        }
        else
        {
            Conditions.Add(condition);
        }

        LastModifiedUtc = DateTime.UtcNow;
    }

    public void AddOrReplaceSpectraEntry(FeatureSpectraEntry spectraEntry)
    {
        ArgumentNullException.ThrowIfNull(spectraEntry);

        int existingIndex = SpectraFiles.FindIndex(p => p.MatchesRawFile(spectraEntry.MassSpecFilePath));
        if (existingIndex >= 0)
        {
            SpectraFiles[existingIndex] = spectraEntry;
        }
        else
        {
            SpectraFiles.Add(spectraEntry);
        }

        LastModifiedUtc = DateTime.UtcNow;
    }

    public bool TryGetSpectraEntry(string rawFilePath, out FeatureSpectraEntry entry)
    {
        entry = SpectraFiles.FirstOrDefault(p => p.MatchesRawFile(rawFilePath));
        return entry is not null;
    }

    public bool TryGetFeatureFile(string rawFilePath, string conditionKey, out string featureFilePath)
    {
        featureFilePath = null;

        if (!TryGetSpectraEntry(rawFilePath, out var spectraEntry))
        {
            return false;
        }

        if (!spectraEntry.TryGetConditionFile(conditionKey, out var conditionFile))
        {
            return false;
        }

        featureFilePath = conditionFile.FeatureFilePath;
        return true;
    }

    public bool TryGetCondition(string conditionKey, out FeatureFileMapCondition condition)
    {
        condition = Conditions.FirstOrDefault(p => p.MatchesKey(conditionKey));
        return condition is not null;
    }

    #endregion

    public SearchFeatureFileMap BuildSearchMap(IEnumerable<string> rawFiles, string conditionKey)
    {
        if (!TryGetCondition(conditionKey, out var condition))
            throw new FeatureMappingException($"Condition '{conditionKey}' not found in map.");

        Dictionary<string, SearchFeatureFileMapEntry?> searchMap = rawFiles.ToDictionary(p => p, _ => (SearchFeatureFileMapEntry)null);

        foreach (var spectraEntry in SpectraFiles.Where(p => searchMap.ContainsKey(p.MassSpecFilePath)))
            if (spectraEntry.TryGetConditionFile(conditionKey, out var conditionFile))
                searchMap[spectraEntry.MassSpecFilePath] = new SearchFeatureFileMapEntry
                {
                    MassSpecFilePath = spectraEntry.MassSpecFilePath,
                    MassSpecFileName = spectraEntry.MassSpecFileName,
                    FeatureFilePath = conditionFile.FeatureFilePath,
                    FeatureFileName = conditionFile.FeatureFileName,
                };

        if (searchMap.Values.Any(p => p is null))
            throw new FeatureMappingException($"The following raw files are missing feature file mappings for condition '{conditionKey}': {string.Join(", ", searchMap.Where(kvp => kvp.Value is null).Select(kvp => kvp.Key))}");

        return new SearchFeatureFileMap
        {
            SourceMapPath = FilePath,
            SelectedConditionKey = condition.ConditionKey,
            SelectedConditionDisplayName = condition.DisplayName,
            Entries = searchMap.Values.Cast<SearchFeatureFileMapEntry>().ToList()
        };
    }
}
