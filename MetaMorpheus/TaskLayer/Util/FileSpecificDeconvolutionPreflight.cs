using MassSpectrometry;
using MzLibUtil;
using Readers;
using System;
using System.Collections.Concurrent;
using System.IO;
using System.Runtime.CompilerServices;
using System.Threading;

namespace TaskLayer.Util;

/// <summary>
/// Performs runtime-only validation for native file-specific FromFile deconvolution parameters.
/// UI validation remains intentionally lightweight; this class forces mzLib's lazy feature load
/// immediately before task processing and memoizes the resulting verdict so a later call cannot
/// contradict the fallback decision (and warning) already made for the same settings and raw file.
/// </summary>
internal static class FileSpecificDeconvolutionPreflight
{
    private static readonly ConditionalWeakTable<FileSpecificParameters, ConcurrentDictionary<string, Lazy<PreflightDecision>>> Decisions = new();

    /// <summary>
    /// Returns the FromFile-or-fallback verdict for this (file-specific settings, raw file) pair,
    /// forcing the lazy feature load at most once. The first call decides; every later call returns
    /// the same verdict without re-attempting the load, so a transient later load cannot diverge
    /// from the fallback already taken and warned about.
    /// </summary>
    internal static PreflightDecision Resolve(
        FileSpecificParameters fileSpecificParams,
        FromFileDeconvolutionParameters fromFileParams,
        string rawFilePath)
    {
        string key = rawFilePath ?? "<unknown raw file>";
        ConcurrentDictionary<string, Lazy<PreflightDecision>> perSettings =
            Decisions.GetOrCreateValue(fileSpecificParams);

        return perSettings.GetOrAdd(
            key,
            _ => new Lazy<PreflightDecision>(
                () => CreateDecision(fromFileParams),
                LazyThreadSafetyMode.ExecutionAndPublication)).Value;
    }

    private static PreflightDecision CreateDecision(FromFileDeconvolutionParameters fromFileParams)
    {
        bool succeeded = TryPreload(fromFileParams, out string failureReason);
        return new PreflightDecision(succeeded, failureReason);
    }

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

    /// <summary>
    /// The memoized outcome of the file-specific preflight. <see cref="TryConsumeWarning"/> hands out
    /// the fallback warning at most once no matter how many callers reuse this decision.
    /// </summary>
    internal sealed class PreflightDecision
    {
        private int _warningConsumed;

        internal PreflightDecision(bool succeeded, string failureReason)
        {
            Succeeded = succeeded;
            FailureReason = failureReason;
        }

        internal bool Succeeded { get; }

        internal string FailureReason { get; }

        /// <summary>Returns true the first time only, so the fallback is logged once per decision.</summary>
        internal bool TryConsumeWarning() => Interlocked.Exchange(ref _warningConsumed, 1) == 0;
    }
}
