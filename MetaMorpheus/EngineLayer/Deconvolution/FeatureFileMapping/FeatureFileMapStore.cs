using Nett;
using System;
using System.IO;

namespace EngineLayer.Deconvolution.FeatureFileMapping;

/// <summary>
/// Master feature-map store service with atomic TOML persistence.
/// Handles loading and saving <see cref="FeatureFileMap"/> to disk,
/// returning a usable empty map when the file does not exist,
/// and throwing <see cref="FeatureMappingException"/> on invalid TOML.
/// </summary>
public static class FeatureFileMapStore
{
    /// <summary>
    /// Normalize a file-system path to its canonical full form.
    /// Resolves relative paths, normalises directory separators to the OS default,
    /// and collapses <c>.</c> / <c>..</c> segments so that equivalent paths have
    /// a single string representation for comparison and storage.
    /// </summary>
    internal static string NormalizePath(string path)
    {
        if (string.IsNullOrEmpty(path))
            return path ?? string.Empty;
        return Path.GetFullPath(path);
    }

    /// <summary>
    /// Load the master feature-map store from the given path.
    /// If <paramref name="storePath" /> is null, <see cref="EngineLayer.GlobalVariables.FeatureMapsFilePath" /> is used.
    /// Returns an empty <see cref="FeatureFileMap" /> when the file does not exist.
    /// </summary>
    /// <param name="storePath">Path to the TOML feature-map file, or null for the default.</param>
    /// <returns>A populated or empty <see cref="FeatureFileMap" />.</returns>
    /// <exception cref="FeatureMappingException">When the file exists but cannot be parsed as a valid feature-map TOML.</exception>
    public static FeatureFileMap Load(string storePath = null)
    {
        storePath = NormalizePath(storePath ?? EngineLayer.GlobalVariables.FeatureMapsFilePath);

        if (!File.Exists(storePath))
        {
            return new FeatureFileMap
            {
                FilePath = storePath,
                LastModifiedUtc = DateTime.UtcNow
            };
        }

        try
        {
            var map = Toml.ReadFile<FeatureFileMap>(storePath);
            map.FilePath = storePath;
            NormalizePathsInMap(map);
            return map;
        }
        catch (Exception ex) when (ex is not FeatureMappingException)
        {
            throw new FeatureMappingException(
                $"Failed to parse feature-map TOML at '{storePath}'. The file may be corrupted or in an invalid format.",
                ex);
        }
    }

    /// <summary>
    /// Apply <see cref="NormalizePath"/> to every path stored inside a loaded map
    /// so that lookups through <see cref="FeatureSpectraEntry.MatchesRawFile"/>
    /// and similar comparisons work with canonical forms.
    /// </summary>
    private static void NormalizePathsInMap(FeatureFileMap map)
    {
        foreach (var spectraEntry in map.SpectraFiles)
        {
            spectraEntry.MassSpecFilePath = NormalizePath(spectraEntry.MassSpecFilePath);
            foreach (var cf in spectraEntry.FeatureFiles)
            {
                cf.FeatureFilePath = NormalizePath(cf.FeatureFilePath);
            }
        }
    }

    /// <summary>
    /// Save the feature-map store atomically via a temp-write-and-replace strategy.
    /// Creates intermediate directories if they do not exist.
    /// </summary>
    /// <param name="map">The feature map to persist.</param>
    /// <param name="storePath">Path to write to, or null for <see cref="EngineLayer.GlobalVariables.FeatureMapsFilePath"/>.</param>
    /// <exception cref="ArgumentNullException">When <paramref name="map"/> is null.</exception>
    public static void Save(FeatureFileMap map, string storePath = null)
    {
        ArgumentNullException.ThrowIfNull(map);

        storePath = NormalizePath(storePath ?? EngineLayer.GlobalVariables.FeatureMapsFilePath);

        map.LastModifiedUtc = DateTime.UtcNow;

        string dir = Path.GetDirectoryName(storePath);
        if (!string.IsNullOrEmpty(dir) && !Directory.Exists(dir))
        {
            Directory.CreateDirectory(dir);
        }

        string tempPath = storePath + ".tmp." + Guid.NewGuid().ToString("N");

        try
        {
            Toml.WriteFile(map, tempPath);

            if (File.Exists(storePath))
            {
                // Atomic replace when destination exists
                File.Replace(tempPath, storePath, null);
            }
            else
            {
                // First save: move temp to destination
                File.Move(tempPath, storePath);
            }
        }
        catch
        {
            // Clean up temp file on failure so no orphaned .tmp.* files linger
            if (File.Exists(tempPath))
            {
                try { File.Delete(tempPath); }
                catch { /* best effort */ }
            }
            throw;
        }
    }
}
