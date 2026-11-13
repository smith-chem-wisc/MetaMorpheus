#nullable enable
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chromatography.RetentionTimeCalibration;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.Chronologer;
using EngineLayer;
using MathNet.Numerics.Statistics;
using MzLibUtil;

namespace TaskLayer.ParallelSearch.Analysis.Collectors;

/// <summary>
/// Compares each confident transient identification's retention time with Chronologer's prediction.
/// Chronologer predicts % acetonitrile, not minutes, so predictions go through a per-file calibration
/// fitted once on the baseline search (<see cref="FitRetentionTimeCalibrations"/>). Every database uses
/// the same calibration, so errors are in minutes and comparable between databases.
/// </summary>
public class RetentionTimeCollector : IMetricCollector
{
    private static readonly RetentionTimePredictor _predictor = new ChronologerRetentionTimePredictor();
    private static readonly ConcurrentDictionary<string, double> _sharedPredictionCache = new(StringComparer.Ordinal);

    public const string PsmMeanAbsoluteRtError = "PsmMeanAbsoluteRtError";
    public const string PsmRtCorrelationCoefficient = "PsmRtCorrelationCoefficient";
    public const string PsmAllRtErrors = "PsmAllRtErrors";
    public const string PeptideMeanAbsoluteRtError = "PeptideMeanAbsoluteRtError";
    public const string PeptideRtCorrelationCoefficient = "PeptideRtCorrelationCoefficient";
    public const string PeptideAllRtErrors = "PeptideAllRtErrors";

    public string CollectorName => "RetentionTime";
    public IEnumerable<string> GetOutputColumns()
    {
        yield return PsmMeanAbsoluteRtError;
        yield return PsmRtCorrelationCoefficient;
        yield return PsmAllRtErrors;
        yield return PeptideMeanAbsoluteRtError;
        yield return PeptideRtCorrelationCoefficient;
        yield return PeptideAllRtErrors;
    }

    /// <summary>
    /// Fits, for each spectra file, Chronologer's predicted % acetonitrile against observed retention time in minutes,
    /// using mzLib's robust LOWESS (<see cref="IrtCalibration.Fit"/>). Anchors are the baseline's confident target PSMs,
    /// one per peptide per file (its best-scoring PSM). A file that cannot be fitted, for example with fewer than 20
    /// anchors, is left out and warned about; its transient PSMs then get no retention time error.
    /// </summary>
    /// <returns>Calibrations keyed by spectra file path.</returns>
    public static Dictionary<string, IrtCalibrationModel> FitRetentionTimeCalibrations(IEnumerable<SpectralMatch> baselinePsms,
        double qValueThreshold, int maxThreads, Action<string> warn)
    {
        var anchorPsms = baselinePsms
            .Where(p => !p.IsDecoy && p.FdrInfo?.QValue <= qValueThreshold && p.FullSequence != null && p.FullFilePath != null)
            .GroupBy(p => (p.FullFilePath, p.FullSequence))
            .Select(g => g.MaxBy(p => p.Score)!)
            .ToList();

        PredictUncachedSequences(anchorPsms, maxThreads);

        var calibrations = new Dictionary<string, IrtCalibrationModel>();
        foreach (var file in anchorPsms.GroupBy(p => p.FullFilePath))
        {
            var anchors = file
                .Select(p => (Rt: new RtMinutes(p.ScanRetentionTime), Irt: new Irt(CachedPrediction(p.FullSequence))))
                .Where(a => double.IsFinite(a.Rt.Value) && double.IsFinite(a.Irt.Value))
                .ToList();
            try
            {
                calibrations[file.Key] = IrtCalibration.Fit(anchors);
            }
            catch (ArgumentException e)
            {
                warn($"Retention time errors are not reported for {Path.GetFileName(file.Key)}: " +
                     $"Chronologer's predictions could not be calibrated to this file's retention times ({e.Message})");
            }
        }

        return calibrations;
    }

    /// <summary>
    /// Predicts, in one batched call, every target sequence among <paramref name="matches"/> that is not cached yet.
    /// One sequence at a time takes the model lock per peptide, so databases searched at once queue on each other.
    /// </summary>
    /// <remarks>
    /// A batched value can differ from another batch's value for the same sequence by float32 rounding (about 1e-5),
    /// so the first value cached wins and every database reads its predictions back from the cache.
    /// </remarks>
    private static void PredictUncachedSequences(IEnumerable<SpectralMatch> matches, int maxThreads)
    {
        var uncached = new Dictionary<string, IRetentionPredictable>(StringComparer.Ordinal);
        foreach (var match in matches)
        {
            foreach (var hypothesis in match.BestMatchingBioPolymersWithSetMods)
            {
                if (hypothesis.IsDecoy || hypothesis.SpecificBioPolymer is not IRetentionPredictable predictable)
                    continue;

                if (!_sharedPredictionCache.ContainsKey(predictable.FullSequence))
                    uncached.TryAdd(predictable.FullSequence, predictable);
            }
        }

        if (uncached.Count == 0)
            return;

        IReadOnlyList<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)> predictions;
        try
        {
            predictions = _predictor.PredictRetentionTimeEquivalents(uncached.Values, maxThreads);
        }
        catch (Exception)
        {
            // A fault in one batched forward pass fails the whole call; the single-peptide path absorbs faults per peptide.
            foreach (var (fullSequence, predictable) in uncached)
                _sharedPredictionCache.TryAdd(fullSequence, _predictor.PredictRetentionTimeEquivalent(predictable, out _) ?? double.NaN);
            return;
        }

        // Keyed by each result's own peptide: the predictor interface does not promise results in input order.
        foreach (var (predictedValue, peptide, _) in predictions)
            _sharedPredictionCache.TryAdd(peptide.FullSequence, predictedValue ?? double.NaN);
    }

    /// <summary>The cached prediction for a sequence, or NaN when it could not be predicted.</summary>
    private static double CachedPrediction(string fullSequence)
    {
        return _sharedPredictionCache.TryGetValue(fullSequence, out double predicted) ? predicted : double.NaN;
    }

    /// <summary>The calibration for the file a match came from, or null when that file has none.</summary>
    private static IrtCalibrationModel? CalibrationFor(TransientDatabaseContext context, SpectralMatch match)
    {
        return match.FullFilePath != null && context.RetentionTimeCalibrations != null
            && context.RetentionTimeCalibrations.TryGetValue(match.FullFilePath, out var calibration)
            ? calibration
            : null;
    }

    /// <summary>A sequence's predicted retention time in minutes on a file's calibration, or NaN when it could not be predicted.</summary>
    private static double PredictedMinutes(IrtCalibrationModel calibration, string fullSequence)
    {
        double predictedAcetonitrile = CachedPrediction(fullSequence);
        return double.IsNaN(predictedAcetonitrile)
            ? double.NaN
            : calibration.ToRtMinutes(new Irt(predictedAcetonitrile)).Value;
    }

    public Dictionary<string, object> CollectData(TransientDatabaseContext context)
    {

        // Psms
        var confidentPsms = context.TransientPsms
            .Where(context.PsmCutoff.Passes)
            .ToList();
        var confidentPeptides = context.TransientPeptides
            .Where(context.PeptideCutoff.Passes)
            .ToList();

        PredictUncachedSequences(confidentPsms.Concat(confidentPeptides), context.CommonParameters.MaxThreadsToUsePerFile);

        List<double> allPsmRtErrors = new();
        List<double> psmObservedRts = new();
        List<double> psmPredictedRts = new();

        foreach (var psm in confidentPsms)
        {
            double observedRt = psm.ScanRetentionTime;
            if (CalibrationFor(context, psm) is not { } calibration)
                continue;

            foreach (var hypothesis in psm.BestMatchingBioPolymersWithSetMods)
            {
                if (hypothesis.IsDecoy)
                    continue;

                double predictedRt = PredictedMinutes(calibration, hypothesis.SpecificBioPolymer.FullSequence);

                if (!double.IsFinite(predictedRt))
                    continue;

                double rtError = observedRt - predictedRt;
                allPsmRtErrors.Add(rtError);
                psmObservedRts.Add(observedRt);
                psmPredictedRts.Add(predictedRt);
                
                break; // Only use first hypothesis for RT prediction
            }
        }

        // Peptides
        List<double> allPeptideRtErrors = new();
        List<double> peptideObservedRts = new();
        List<double> peptidePredictedRts = new();

        foreach (var peptide in confidentPeptides)
        {
            double observedRt = peptide.ScanRetentionTime;
            if (CalibrationFor(context, peptide) is not { } calibration)
                continue;

            foreach (var hypothesis in peptide.BestMatchingBioPolymersWithSetMods)
            {
                if (hypothesis.IsDecoy)
                    continue;

                double predictedRt = PredictedMinutes(calibration, hypothesis.SpecificBioPolymer.FullSequence);

                if (!double.IsFinite(predictedRt))
                    continue;

                double rtError = observedRt - predictedRt;
                allPeptideRtErrors.Add(rtError);
                peptideObservedRts.Add(observedRt);
                peptidePredictedRts.Add(predictedRt);
                
                break; // Only use first hypothesis for RT prediction
            }
        }

        // Calculate statistics - USE NaN FOR INSUFFICIENT DATA
        double psmMeanAbsoluteError = allPsmRtErrors.Any()
            ? allPsmRtErrors.Select(Math.Abs).Mean()
            : double.NaN; // Changed from 0 to NaN

        double psmCorrelation = psmObservedRts.Count > 1
            ? Correlation.Pearson(psmObservedRts, psmPredictedRts)
            : double.NaN; // Changed from 0 to NaN

        double peptideMeanAbsoluteError = allPeptideRtErrors.Any()
            ? allPeptideRtErrors.Select(Math.Abs).Mean()
            : double.NaN; // Changed from 0 to NaN

        double peptideCorrelation = peptideObservedRts.Count > 1
            ? Correlation.Pearson(peptideObservedRts, peptidePredictedRts)
            : double.NaN; // Changed from 0 to NaN

        return new Dictionary<string, object>
        {
            [PsmMeanAbsoluteRtError] = psmMeanAbsoluteError,
            [PsmRtCorrelationCoefficient] = psmCorrelation,
            [PsmAllRtErrors] = allPsmRtErrors.ToArray(),
            [PeptideMeanAbsoluteRtError] = peptideMeanAbsoluteError,
            [PeptideRtCorrelationCoefficient] = peptideCorrelation,
            [PeptideAllRtErrors] = allPeptideRtErrors.ToArray()
        };
    }

    public bool CanCollectData(TransientDatabaseContext context)
    {
        return context.AllPsms != null
               && context.TransientPsms != null
               && context.TransientPsms.All(p => p is PeptideSpectralMatch)
               && context.AllPeptides != null
               && context.TransientPeptides != null
               && context.TransientPeptides.All(p => p is PeptideSpectralMatch);
    }
}
