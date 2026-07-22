using Nett;
using System;
using System.Collections.Generic;
using System.Linq;
using EngineLayer.Deconvolution.FeatureFileMapping;

namespace EngineLayer.Deconvolution;

/// <summary>
/// Feature file mapping used at time of search. This is a simplified version of the FeatureFileMap used for deconvolution, containing only the necessary information for searching.
/// </summary>
public class SearchFeatureFileMap : IEquatable<SearchFeatureFileMap>
{
    public string SourceMapPath { get; set; } = string.Empty;
    public string SelectedConditionKey { get; set; } = string.Empty;
    public string SelectedConditionDisplayName { get; set; } = string.Empty;
    public List<SearchFeatureFileMapEntry> Entries { get; set; } = new();

    public bool TryGetFeaturePathForMassSpecFile(string massSpecFilePath, out string featureFilePath)
    {
        featureFilePath = string.Empty;
        var entry = Entries.Find(e => string.Equals(e.MassSpecFilePath, massSpecFilePath, StringComparison.OrdinalIgnoreCase));
        if (entry != null)
        {
            featureFilePath = entry.FeatureFilePath;
            return true;
        }
        return false;
    }

    public SearchFeatureFileMap Clone()
    {
        return new SearchFeatureFileMap
        {
            SourceMapPath = this.SourceMapPath,
            SelectedConditionKey = this.SelectedConditionKey,
            SelectedConditionDisplayName = this.SelectedConditionDisplayName,
            Entries = new List<SearchFeatureFileMapEntry>(this.Entries.Select(p => p.Clone()))
        };
    }

    /// <summary>
    /// Returns true when this map has no entries and would fail at materialization time.
    /// </summary>
    [TomlIgnore]
    public bool IsEmpty => Entries == null || Entries.Count == 0;

    /// <summary>
    /// Throws <see cref="FeatureMappingException"/> if the map is empty (no entries).
    /// Call this before materialization to fail early with a clear message.
    /// </summary>
    public void ValidateNotEmpty()
    {
        if (IsEmpty)
        {
            throw new FeatureMappingException(
                "Search feature file map is empty: no feature file entries are embedded in the task. " +
                "Deconvolution cannot proceed without at least one feature-file-to-spectra-file mapping.");
        }
    }

    public bool Equals(SearchFeatureFileMap other)
    {
        if (other is null)
        {
            return false;
        }

        if (ReferenceEquals(this, other))
        {
            return true;
        }

        // SourceMapPath is audit-only metadata; it does not affect execution equality.
        return string.Equals(SelectedConditionKey, other.SelectedConditionKey, StringComparison.Ordinal)
            && string.Equals(SelectedConditionDisplayName, other.SelectedConditionDisplayName, StringComparison.Ordinal)
            && Entries.SequenceEqual(other.Entries);
    }

    public override bool Equals(object obj) => Equals(obj as SearchFeatureFileMap);

    public override int GetHashCode()
    {
        var hash = new HashCode();
        // SourceMapPath is intentionally excluded: it is audit-only metadata.
        hash.Add(SelectedConditionKey, StringComparer.Ordinal);
        hash.Add(SelectedConditionDisplayName, StringComparer.Ordinal);
        foreach (var entry in Entries)
        {
            hash.Add(entry);
        }
        return hash.ToHashCode();
    }
}

/// <summary>
/// Represents a single entry in the SearchFeatureFileMap, associating a mass spectrometry file with its corresponding feature file for a specific condition.
/// </summary>
public class SearchFeatureFileMapEntry : IEquatable<SearchFeatureFileMapEntry>
{
    public string MassSpecFilePath { get; set; } = string.Empty;
    public string MassSpecFileName { get; set; } = string.Empty;
    public string FeatureFilePath { get; set; } = string.Empty;
    public string FeatureFileName { get; set; } = string.Empty;

    public SearchFeatureFileMapEntry Clone()
    {
        return new SearchFeatureFileMapEntry
        {
            MassSpecFilePath = MassSpecFilePath,
            MassSpecFileName = MassSpecFileName,
            FeatureFilePath = FeatureFilePath,
            FeatureFileName = FeatureFileName,
        };
    }

    public bool Equals(SearchFeatureFileMapEntry other)
    {
        if (other is null)
        {
            return false;
        }

        if (ReferenceEquals(this, other))
        {
            return true;
        }

        return string.Equals(MassSpecFilePath, other.MassSpecFilePath, StringComparison.Ordinal)
            && string.Equals(MassSpecFileName, other.MassSpecFileName, StringComparison.Ordinal)
            && string.Equals(FeatureFilePath, other.FeatureFilePath, StringComparison.Ordinal)
            && string.Equals(FeatureFileName, other.FeatureFileName, StringComparison.Ordinal);
    }

    public override bool Equals(object obj) => Equals(obj as SearchFeatureFileMapEntry);

    public override int GetHashCode()
        => HashCode.Combine(MassSpecFilePath, MassSpecFileName, FeatureFilePath, FeatureFileName);

    public override string ToString() => $"{MassSpecFilePath},{FeatureFilePath}";
}
