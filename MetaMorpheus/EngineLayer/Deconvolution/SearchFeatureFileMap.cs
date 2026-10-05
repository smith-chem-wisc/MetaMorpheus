using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Deconvolution;

/// <summary>
/// Explicit raw-file to MS1 feature-file mappings embedded in task parameters.
/// </summary>
public class SearchFeatureFileMap : IEquatable<SearchFeatureFileMap>
{
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
            Entries = new List<SearchFeatureFileMapEntry>(this.Entries.Select(p => p.Clone()))
        };
    }

    /// <summary>
    /// Returns true when no raw-file mappings are configured.
    /// </summary>
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

        return Entries.SequenceEqual(other.Entries);
    }

    public override bool Equals(object obj) => Equals(obj as SearchFeatureFileMap);

    public override int GetHashCode()
    {
        var hash = new HashCode();
        foreach (var entry in Entries)
        {
            hash.Add(entry);
        }
        return hash.ToHashCode();
    }
}

/// <summary>
/// Associates one mass spectrometry file with its corresponding MS1 feature file.
/// </summary>
public class SearchFeatureFileMapEntry : IEquatable<SearchFeatureFileMapEntry>
{
    public string MassSpecFilePath { get; set; } = string.Empty;
    public string FeatureFilePath { get; set; } = string.Empty;

    public SearchFeatureFileMapEntry Clone()
    {
        return new SearchFeatureFileMapEntry
        {
            MassSpecFilePath = MassSpecFilePath,
            FeatureFilePath = FeatureFilePath,
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
            && string.Equals(FeatureFilePath, other.FeatureFilePath, StringComparison.Ordinal);
    }

    public override bool Equals(object obj) => Equals(obj as SearchFeatureFileMapEntry);

    public override int GetHashCode()
        => HashCode.Combine(MassSpecFilePath, FeatureFilePath);

    public override string ToString() => $"{MassSpecFilePath},{FeatureFilePath}";
}
