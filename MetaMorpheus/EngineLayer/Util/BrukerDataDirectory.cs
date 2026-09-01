using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MzLibUtil;
using Readers;

namespace EngineLayer.Util
{
    /// <summary>
    /// Single source of truth for "is this path Bruker data that MetaMorpheus can read, and where should the reader be
    /// pointed?". Bruker acquisitions are directories ending in ".d"; the flavour is decided by the files inside:
    ///
    ///     qTOF        analysis.baf                      -> mzLib BrukerFileReader
    ///     timsTOF     analysis.tdf + analysis.tdf_bin   -> mzLib TimsTofFileReader (TDF)
    ///     TIMS off    analysis.tsf + analysis.tsf_bin   -> mzLib TimsTofFileReader (TSF)
    ///
    /// Classification is delegated to mzLib's ParseFileType, which is the same call MsDataFileReader.GetDataFile uses to
    /// choose a reader. Re-implementing the check here would let MetaMorpheus's front door drift from what mzLib can
    /// actually open, which is exactly how ".tsf" support went missing after mzLib gained it.
    /// </summary>
    public static class BrukerDataDirectory
    {
        public const string DotD = ".d";

        /// <summary>
        /// Extensions of the individual files that live inside a ".d" folder. Users routinely hand us one of these
        /// instead of the folder itself, so every entry point has to redirect them to the parent ".d".
        /// </summary>
        public static readonly IReadOnlyList<string> InnerFileExtensions =
            new[] { ".baf", ".tdf", ".tdf_bin", ".tsf", ".tsf_bin" };

        /// <summary>
        /// True if the extension belongs to a file found inside a Bruker ".d" folder.
        /// </summary>
        public static bool IsInnerFileExtension(string extension)
        {
            return extension != null && InnerFileExtensions.Contains(extension.ToLowerInvariant());
        }

        /// <summary>
        /// True if the path names a ".d" folder. Says nothing about whether the folder holds readable data.
        /// </summary>
        public static bool IsDotDPath(string path)
        {
            return path != null && path.EndsWith(DotD, StringComparison.OrdinalIgnoreCase);
        }

        /// <summary>
        /// True if the path is a ".d" folder holding Bruker data that mzLib can open. Never throws: a folder that is not
        /// Bruker data is an expected answer here, not a fault.
        /// </summary>
        public static bool IsValid(string directoryPath)
        {
            if (!IsDotDPath(directoryPath) || !Directory.Exists(directoryPath))
            {
                return false;
            }

            SupportedFileType fileType;
            try
            {
                fileType = directoryPath.ParseFileType();
            }
            catch (Exception e) when (e is MzLibException || e is IOException || e is UnauthorizedAccessException)
            {
                return false; // not Bruker data that mzLib recognizes
            }

            switch (fileType)
            {
                case SupportedFileType.BrukerD:
                    return true;

                // ParseFileType classifies timsTOF from the SQLite file alone, but TimsTofFileReader also needs the
                // binary sidecar - InitiateDynamicConnection throws FileNotFoundException without it. Requiring it here
                // keeps "MetaMorpheus says this is readable" and "mzLib can actually open it" in agreement.
                case SupportedFileType.BrukerTimsTof:
                    return HasTimsBinarySidecar(directoryPath);

                default:
                    return false;
            }
        }

        private static bool HasTimsBinarySidecar(string directoryPath)
        {
            return (File.Exists(Path.Combine(directoryPath, "analysis.tdf"))
                        && File.Exists(Path.Combine(directoryPath, "analysis.tdf_bin")))
                || (File.Exists(Path.Combine(directoryPath, "analysis.tsf"))
                        && File.Exists(Path.Combine(directoryPath, "analysis.tsf_bin")));
        }

        /// <summary>
        /// If the path is a file inside a readable Bruker ".d" folder, hands back that folder. mzLib's Bruker readers are
        /// pointed at the directory, never at the inner file, so every place that accepts a spectra path has to make this
        /// substitution.
        /// </summary>
        /// <returns>true if <paramref name="dotDFolder"/> was set; false if the path is not a redirectable inner file.</returns>
        public static bool TryGetParentDotDFolder(string path, out string dotDFolder)
        {
            dotDFolder = null;

            if (string.IsNullOrEmpty(path) || !IsInnerFileExtension(Path.GetExtension(path)))
            {
                return false;
            }

            string parentDirectory = Path.GetDirectoryName(path);
            if (!IsValid(parentDirectory))
            {
                return false;
            }

            dotDFolder = parentDirectory;
            return true;
        }
    }
}
