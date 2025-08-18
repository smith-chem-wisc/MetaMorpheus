using System;
using System.IO;
using System.Linq;

namespace EngineLayer.Util;

public static class PathSafety
{
    // Conservative defaults for broad Windows tooling compatibility.
    // You can raise these if your environment supports long paths (\\?\).
    public const int DefaultMaxPath = 260;
    public const int DefaultMaxFileName = 255;

    /// <summary>
    /// Validates and sanitizes a path for an output file. Ensures the result:
    /// - Ends with <paramref name="requiredEnding"/> (case-insensitive; added if missing)
    /// - Has a filename with invalid characters replaced by '_'
    /// - Avoids Windows reserved device names by prefixing an underscore
    /// - Trims the base filename (before the required ending) to fit within path and filename limits
    /// Throws if the directory portion alone exceeds <paramref name="maxPath"/> and cannot be fixed.
    /// </summary>
    /// <param name="pathToValidate">Proposed output path (may be relative or absolute).</param>
    /// <param name="requiredEnding">Required trailing suffix (e.g., ".mzML"). Case-insensitive.</param>
    /// <param name="maxPath">Max total path length to enforce (default 260).</param>
    /// <param name="maxFileName">Max filename length to enforce (default 255).</param>
    /// <returns>Sanitized, safe path.</returns>
    /// <exception cref="ArgumentException">On null/empty inputs or invalid ending.</exception>
    /// <exception cref="PathTooLongException">If the directory part alone makes the path exceed limits.</exception>
    public static string MakeSafeOutputPath(
        string pathToValidate,
        string requiredEnding,
        int maxPath = DefaultMaxPath,
        int maxFileName = DefaultMaxFileName)
    {
        if (string.IsNullOrWhiteSpace(pathToValidate))
            throw new ArgumentException("Path is null or whitespace.", nameof(pathToValidate));
        if (string.IsNullOrWhiteSpace(requiredEnding))
            throw new ArgumentException("Required ending is null or whitespace.", nameof(requiredEnding));

        // Normalize required ending to begin with '.'
        if (!requiredEnding.StartsWith(".", StringComparison.Ordinal))
            requiredEnding = "." + requiredEnding;
        // Use a consistent case for comparison, but we’ll preserve the supplied ending’s case when appending.
        var reqEndingCmp = requiredEnding.ToLowerInvariant();

        // Split directory + filename
        string directory = Path.GetDirectoryName(pathToValidate) ?? string.Empty;
        string fileName = Path.GetFileName(pathToValidate);

        // If no filename provided (path ends with directory), create a placeholder filename
        if (string.IsNullOrWhiteSpace(fileName))
            fileName = "output" + requiredEnding;

        // Ensure the filename ends with required ending (case-insensitive). If not, append it.
        if (!fileName.ToLowerInvariant().EndsWith(reqEndingCmp))
            fileName += requiredEnding;

        // Sanitize filename: replace invalid characters
        var invalid = Path.GetInvalidFileNameChars();
        var sanitizedName = new string(fileName.Select(c => invalid.Contains(c) ? '_' : c).ToArray());

        // Ensure base name not reserved device (CON, PRN, AUX, NUL, COM1..9, LPT1..9)
        sanitizedName = AvoidReservedDeviceNames(sanitizedName);

        // Enforce filename length limit by trimming the base (before required ending), preserving the ending
        // Identify the ending position case-insensitively
        var lower = sanitizedName.ToLowerInvariant();
        int endingIndex = lower.LastIndexOf(reqEndingCmp, StringComparison.Ordinal);
        string basePart = sanitizedName.Substring(0, endingIndex);
        string endingPart = sanitizedName.Substring(endingIndex); // includes the dot and suffix

        // If filename too long, trim the base to fit maxFileName
        if (sanitizedName.Length > maxFileName)
        {
            int maxBaseLen = Math.Max(1, maxFileName - endingPart.Length);
            if (basePart.Length > maxBaseLen)
                basePart = basePart.Substring(0, maxBaseLen);
            sanitizedName = basePart + endingPart;
        }

        // Recombine with directory
        string combined = CombineDirectoryAndFile(directory, sanitizedName);

        // Enforce overall path limit by trimming base further if necessary
        if (combined.Length > maxPath)
        {
            // How many chars can filename have given directory length + path separator
            string dirWithSep = EnsureDirWithSeparator(directory);
            int maxFileLenGivenDir = maxPath - dirWithSep.Length;
            if (maxFileLenGivenDir <= 0)
                throw new PathTooLongException("Directory portion exceeds maximum path length limit.");

            // Trim base to fit the remaining allowed filename length
            int allowedBase = Math.Max(1, Math.Min(basePart.Length, maxFileLenGivenDir - endingPart.Length));
            if (allowedBase <= 0)
                throw new PathTooLongException("Cannot construct a valid filename within the path length limit.");

            basePart = basePart.Substring(0, allowedBase);
            sanitizedName = basePart + endingPart;
            combined = dirWithSep + sanitizedName;

            // Final defensive check
            if (combined.Length > maxPath)
                throw new PathTooLongException("Resulting path exceeds maximum path length after trimming.");
        }

        return combined;
    }

    private static string EnsureDirWithSeparator(string directory)
    {
        if (string.IsNullOrEmpty(directory))
            return string.Empty;
        char sep = Path.DirectorySeparatorChar;
        return directory.EndsWith(sep.ToString(), StringComparison.Ordinal) ? directory : directory + sep;
    }

    private static string CombineDirectoryAndFile(string directory, string fileName)
    {
        if (string.IsNullOrEmpty(directory))
            return fileName;
        return Path.Combine(directory, fileName);
    }

    private static string AvoidReservedDeviceNames(string fileName)
    {
        // Extract base (before extension) for checking reserved names
        string nameNoExt = Path.GetFileNameWithoutExtension(fileName);
        string ext = Path.GetExtension(fileName);

        // Windows device names (case-insensitive)
        // CON, PRN, AUX, NUL, COM1..COM9, LPT1..LPT9
        var reserved = new[]
        {
            "con", "prn", "aux", "nul",
            "com1","com2","com3","com4","com5","com6","com7","com8","com9",
            "lpt1","lpt2","lpt3","lpt4","lpt5","lpt6","lpt7","lpt8","lpt9"
        };

        if (reserved.Contains(nameNoExt.ToLowerInvariant()))
            nameNoExt = "_" + nameNoExt;

        return nameNoExt + ext;
    }
}
