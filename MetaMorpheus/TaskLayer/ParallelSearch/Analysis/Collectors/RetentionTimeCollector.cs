using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.Chronologer;
using EngineLayer;
using MathNet.Numerics.Statistics;

namespace TaskLayer.ParallelSearch.Analysis.Collectors;

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

    private static double GetOrPredictRt(string fullSequence, IRetentionPredictable predictable)
    {
        return _sharedPredictionCache.GetOrAdd(fullSequence, _ =>
        {
            double? predicted = _predictor.PredictRetentionTimeEquivalent(predictable, out var failureReason);
            return predicted.GetValueOrDefault(-1);
        });
    }

    public Dictionary<string, object> CollectData(TransientDatabaseContext context)
    {
        double qValueThreshold = Math.Min(context.CommonParameters.QValueThreshold, context.CommonParameters.PepQValueThreshold);

        // Psms
        var confidentPsms = context.TransientPsms
            .Where(p => p.FdrInfo.QValue <= qValueThreshold)
            .ToList();

        List<double> allPsmRtErrors = new();
        List<double> psmObservedRts = new();
        List<double> psmPredictedRts = new();

        foreach (var psm in confidentPsms)
        {
            double observedRt = psm.ScanRetentionTime;
            
            foreach (var hypothesis in psm.BestMatchingBioPolymersWithSetMods)
            {
                if (hypothesis.IsDecoy)
                    continue;

                string key = hypothesis.SpecificBioPolymer.FullSequence;
                double predictedRt = GetOrPredictRt(key, hypothesis.SpecificBioPolymer as IRetentionPredictable);

                if (predictedRt <= 0)
                    continue;

                double rtError = observedRt - predictedRt;
                allPsmRtErrors.Add(rtError);
                psmObservedRts.Add(observedRt);
                psmPredictedRts.Add(predictedRt);
                
                break; // Only use first hypothesis for RT prediction
            }
        }

        // Peptides
        var confidentPeptides = context.TransientPeptides
            .Where(p => p.FdrInfo.QValue <= qValueThreshold)
            .ToList();

        List<double> allPeptideRtErrors = new();
        List<double> peptideObservedRts = new();
        List<double> peptidePredictedRts = new();

        foreach (var peptide in confidentPeptides)
        {
            double observedRt = peptide.ScanRetentionTime;
            
            foreach (var hypothesis in peptide.BestMatchingBioPolymersWithSetMods)
            {
                if (hypothesis.IsDecoy)
                    continue;

                string key = hypothesis.SpecificBioPolymer.FullSequence;
                double predictedRt = GetOrPredictRt(key, hypothesis.SpecificBioPolymer as IRetentionPredictable);

                if (predictedRt <= 0)
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
