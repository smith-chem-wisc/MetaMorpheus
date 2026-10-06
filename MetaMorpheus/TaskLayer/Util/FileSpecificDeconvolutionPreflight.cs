using MassSpectrometry;
using MzLibUtil;
using Readers;
using System;
using System.Collections.Concurrent;
using System.IO;
using System.Runtime.CompilerServices;

namespace TaskLayer.Util;

/// <summary>
/// Performs runtime-only validation for native file-specific FromFile deconvolution parameters.
/// UI validation remains intentionally lightweight; this class forces mzLib's lazy feature load
/// immediately before task processing and keeps warning deduplication scoped to parsed settings.
/// </summary>
internal static class FileSpecificDeconvolutionPreflight
{
    private static readonly ConditionalWeakTable<FileSpecificParameters, ConcurrentDictionary<string, byte>> WarnedFallbackRawPaths = new();

    internal static bool ShouldWarnFor(FileSpecificParameters fileSpecificParams, string rawFilePath) =>
        WarnedFallbackRawPaths.GetOrCreateValue(fileSpecificParams).TryAdd(rawFilePath ?? "<unknown raw file>", 0);

    /// <summary>
    /// Forces the lazy feature load so malformed or empty sources are diagnosed before scan processing.
    /// Only expected reader, IO, and parsing failures are converted into a failure reason.
    /// </summary>
    internal static bool TryPreload(
        FromFileDeconvolutionParameters fromFileParams,
        out string failureReason)
    {
        failureReason = null;

        if (string.IsNullOrWhiteSpace(fromFileParams.FilePath))
        {
            failureReason = "no feature file path is set";
            return false;
        }

        try
        {
            if (fromFileParams.Features.Count == 0)
            {
                failureReason = "the feature file contains no MS1 features";
                return false;
            }

            return true;
        }
        catch (FileNotFoundException)
        {
            failureReason = "the feature file does not exist";
        }
        catch (DirectoryNotFoundException)
        {
            failureReason = "the feature file's directory does not exist";
        }
        catch (UnauthorizedAccessException)
        {
            failureReason = "access to the feature file was denied";
        }
        catch (IOException e)
        {
            failureReason = $"the feature file could not be read ({e.Message})";
        }
        catch (CsvHelper.CsvHelperException e)
        {
            failureReason = $"the feature file could not be parsed ({e.Message})";
        }
        catch (MzLibException e)
        {
            failureReason = e.Message;
        }

        return false;
    }
}
