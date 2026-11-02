#nullable enable
using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;

namespace EngineLayer.DatabaseLoading;

/// <summary>
/// Reads the first deflines of a FASTA database, compressed or not, so a header setting can be checked
/// against the file before it is loaded.
/// </summary>
public static class FastaHeaderSampler
{
    public const int DefaultHeaderCount = 5;

    /// <summary>
    /// Whether <see cref="DatabaseLoadingEngine"/> would load this path as FASTA, by the same extension rule.
    /// </summary>
    public static bool IsFasta(string path) => InnerExtension(path) is ".fasta" or ".fa";

    /// <summary>
    /// Up to <paramref name="count"/> lines starting with '>', or null when the file cannot be read.
    /// </summary>
    public static List<string>? ReadFirstHeaders(string path, int count = DefaultHeaderCount)
    {
        try
        {
            using Stream file = File.OpenRead(path);
            using Stream content = IsCompressed(path) ? new GZipStream(file, CompressionMode.Decompress) : file;
            using var reader = new StreamReader(content);

            var headers = new List<string>();
            while (headers.Count < count && reader.ReadLine() is { } line)
            {
                if (line.StartsWith('>'))
                    headers.Add(line);
            }
            return headers;
        }
        catch (Exception e) when (e is IOException or InvalidDataException or UnauthorizedAccessException)
        {
            return null;
        }
    }

    // Same rule as the loader: anything ending in "gz" (.gz, .bgz, .tgz) is decompressed.
    private static bool IsCompressed(string path) => Path.GetExtension(path).ToLowerInvariant().EndsWith("gz");

    private static string InnerExtension(string path) => IsCompressed(path)
        ? Path.GetExtension(Path.GetFileNameWithoutExtension(path)).ToLowerInvariant()
        : Path.GetExtension(path).ToLowerInvariant();
}
