using System;
using System.IO;
using System.Linq;

namespace EngineLayer.Util
{
    public static class PathSafety
    {
        // Conservative defaults for broad Windows tooling compatibility.
        // Increase if your environment supports long paths (\\?\).
        public const int DefaultMaxPath = 260;
        public const int DefaultMaxFileName = 255;

        /// <summary>
        /// Validates and sanitizes a path for an output file. Ensures the result:
        /// - Ends with <paramref name="requiredEnding"/> (case-insensitive; appended if missing). The ending may be a complex suffix,
        ///   e.g., "-{identifier}.ext", not just a simple extension.
        /// - Replaces invalid filename characters in the base portion with '_', preserving the required ending verbatim.
        /// - Avoids Windows reserved device names by prefixing an underscore (applied to the name-without-extension).
        /// - Trims only the base portion (before the required ending) to respect filename and total path limits.
        /// Throws if the directory portion alone exceeds <paramref name="maxPath"/> or if a valid filename cannot be constructed.
        /// </summary>
        /// <param name="pathToValidate">Proposed output path (relative or absolute).</param>
        /// <param name="requiredEnding">
        /// Required trailing suffix to enforce (verbatim). Examples: ".mzML", "-run42.mzML", "_v1.json".
        /// Matching is case-insensitive, but the supplied suffix text is preserved on output.
        /// </param>
        /// <param name="maxPath">Maximum total path length (default 260).</param>
        /// <param name="maxFileName">Maximum filename length (default 255).</param>
        /// <returns>Sanitized, safe path.</returns>
        /// <exception cref="ArgumentException">Null/empty inputs.</exception>
        /// <exception cref="PathTooLongException">Directory alone exceeds limit, or no valid filename can be constructed within limits.</exception>
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

            // Split directory and filename
            string directory = Path.GetDirectoryName(pathToValidate) ?? string.Empty;
            string fileName = Path.GetFileName(pathToValidate);

            // If no filename provided, synthesize a base
            if (string.IsNullOrWhiteSpace(fileName))
                fileName = "output";

            // Ensure the filename ends with the requiredEnding (case-insensitive), without duplicating.
            // We will preserve the caller's requiredEnding verbatim.
            bool hasEnding = fileName.EndsWith(requiredEnding, StringComparison.OrdinalIgnoreCase);
            string basePart;
            if (hasEnding)
            {
                // Split into base + ending using the requiredEnding length from the end.
                basePart = fileName.Substring(0, fileName.Length - requiredEnding.Length);
            }
            else
            {
                basePart = fileName;
            }

            // Sanitize only the base portion (so we keep the requiredEnding exactly as provided).
            var invalid = Path.GetInvalidFileNameChars();
            basePart = new string(basePart.Select(c => invalid.Contains(c) ? '_' : c).ToArray());

            // Reassemble filename, appending ending only if it wasn't already present.
            string endingPart = requiredEnding;
            string combinedFileName = basePart + (hasEnding ? fileName[^requiredEnding.Length..] : endingPart);

            // Avoid reserved device names (Windows) by checking the name without extension
            combinedFileName = AvoidReservedDeviceNames(combinedFileName);

            // Enforce filename length by trimming base (not the required ending)
            if (combinedFileName.Length > maxFileName)
            {
                int maxBaseLen = Math.Max(1, maxFileName - endingPart.Length);
                // Recompute basePart against the current combinedFileName to be safe:
                string currentBase = combinedFileName.Substring(0, combinedFileName.Length - endingPart.Length);
                if (currentBase.Length > maxBaseLen)
                    currentBase = currentBase.Substring(0, maxBaseLen);

                combinedFileName = currentBase + endingPart;
            }

            // Combine with directory
            string result = CombineDirectoryAndFile(directory, combinedFileName);

            // Enforce overall path limit by trimming base further if necessary
            if (result.Length > maxPath)
            {
                string dirWithSep = EnsureDirWithSeparator(directory);
                int maxFileLenGivenDir = maxPath - dirWithSep.Length;
                if (maxFileLenGivenDir <= 0)
                    throw new PathTooLongException("Directory portion exceeds maximum path length limit.");

                int allowedBase = Math.Max(1, maxFileLenGivenDir - endingPart.Length);
                if (allowedBase <= 0)
                    throw new PathTooLongException("Cannot construct a valid filename within the path length limit.");

                // Trim base accordingly
                string trimmedBase = combinedFileName.Substring(0, Math.Max(1, allowedBase));
                if (trimmedBase.Length > allowedBase)
                    trimmedBase = trimmedBase.Substring(0, allowedBase);

                combinedFileName = trimmedBase + endingPart;
                result = dirWithSep + combinedFileName;

                if (result.Length > maxPath)
                    throw new PathTooLongException("Resulting path exceeds maximum path length after trimming.");
            }

            return result;
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
            string nameNoExt = Path.GetFileNameWithoutExtension(fileName);
            string ext = Path.GetExtension(fileName);

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
}
