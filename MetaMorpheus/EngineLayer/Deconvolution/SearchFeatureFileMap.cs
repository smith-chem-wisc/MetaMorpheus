using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Nett;

namespace EngineLayer.Deconvolution;

/// <summary>
/// Explicit raw-file to MS1 feature-file mappings embedded in task parameters.
/// </summary>
public class SearchFeatureFileMap : Dictionary<string, string>, IEquatable<SearchFeatureFileMap>
{
    public SearchFeatureFileMap() : base(StringComparer.OrdinalIgnoreCase)
    {
    }

    public SearchFeatureFileMap(IEnumerable<SearchFeatureFileMapEntry> entries) : this()
    {
        foreach (var entry in entries)
            Add(entry.MassSpecFilePath, entry.FeatureFilePath);
    }

    public static SearchFeatureFileMap FromString(string tomlTableContent)
    {
        TomlTable table = Toml.ReadString<TomlTable>(tomlTableContent);
        var map = new SearchFeatureFileMap();
        foreach (string key in table.Keys)
            map.Add(key, table.Get<string>(key));
        return map;
    }

    public override string ToString()
        => string.Join(Environment.NewLine, this.Select(entry =>
            $"{QuoteTomlString(entry.Key)} = {QuoteTomlString(entry.Value)}"));

    private static string QuoteTomlString(string value)
        => $"\"{value.Replace("\\", "\\\\", StringComparison.Ordinal).Replace("\"", "\\\"", StringComparison.Ordinal).Replace("\t", "\\t", StringComparison.Ordinal).Replace("\r", "\\r", StringComparison.Ordinal).Replace("\n", "\\n", StringComparison.Ordinal)}\"";

    public bool TryGetFeaturePathForMassSpecFile(string massSpecFilePath, out string featureFilePath)
    {
        featureFilePath = string.Empty;
        if (TryGetValue(massSpecFilePath, out var exactPath))
        {
            featureFilePath = ResolveFeaturePath(massSpecFilePath, exactPath);
            return true;
        }

        string fileName = Path.GetFileName(massSpecFilePath);
        string fileNameWithoutExtension = Path.GetFileNameWithoutExtension(massSpecFilePath);
        if (TryGetValue(fileName, out var featureFile)
            || TryGetValue(fileNameWithoutExtension, out featureFile))
        {
            featureFilePath = ResolveFeaturePath(massSpecFilePath, featureFile);
            return true;
        }
        return false;
    }

    private static string ResolveFeaturePath(string massSpecFilePath, string featureFilePath)
        => Path.IsPathRooted(featureFilePath)
            ? featureFilePath
            : Path.GetFullPath(Path.Combine(Path.GetDirectoryName(Path.GetFullPath(massSpecFilePath)), featureFilePath));

    public SearchFeatureFileMap Clone()
    {
        return new SearchFeatureFileMap(this.Select(p => new SearchFeatureFileMapEntry(p.Key, p.Value)));
    }

    /// <summary>
    /// Returns true when no raw-file mappings are configured.
    /// </summary>
    public bool IsEmpty => Count == 0;

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

        return Count == other.Count
            && this.All(entry => other.TryGetValue(entry.Key, out var otherFeaturePath)
                && string.Equals(entry.Value, otherFeaturePath, StringComparison.Ordinal));
    }

    public override bool Equals(object obj) => Equals(obj as SearchFeatureFileMap);

    public override int GetHashCode()
    {
        var hash = new HashCode();
        foreach (var entry in this.OrderBy(entry => entry.Key, StringComparer.OrdinalIgnoreCase))
        {
            hash.Add(entry.Key, StringComparer.OrdinalIgnoreCase);
            hash.Add(entry.Value, StringComparer.Ordinal);
        }
        return hash.ToHashCode();
    }
}

/// <summary>
/// Associates one mass spectrometry file with its corresponding MS1 feature file.
/// </summary>
public class SearchFeatureFileMapEntry : IEquatable<SearchFeatureFileMapEntry>
{
    public string MassSpecFilePath { get; set; }
    public string FeatureFilePath { get; set; }

    public SearchFeatureFileMapEntry(string massSpecFilePath, string featureFilePath)
    {
        if (Path.IsPathRooted(massSpecFilePath))
        {
            string fullMassSpecFilePath = Path.GetFullPath(massSpecFilePath);
            MassSpecFilePath = Path.GetFileName(fullMassSpecFilePath);
            FeatureFilePath = Path.IsPathRooted(featureFilePath)
                ? NormalizeTomlPath(Path.GetRelativePath(Path.GetDirectoryName(fullMassSpecFilePath), Path.GetFullPath(featureFilePath)))
                : NormalizeTomlPath(featureFilePath);
        }
        else
        {
            MassSpecFilePath = massSpecFilePath;
            FeatureFilePath = NormalizeTomlPath(featureFilePath);
        }
    }

    private static string NormalizeTomlPath(string path)
        => path.Replace(Path.DirectorySeparatorChar, '/').Replace(Path.AltDirectorySeparatorChar, '/');

    public SearchFeatureFileMapEntry Clone() => new SearchFeatureFileMapEntry(massSpecFilePath: MassSpecFilePath, featureFilePath: FeatureFilePath);

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
