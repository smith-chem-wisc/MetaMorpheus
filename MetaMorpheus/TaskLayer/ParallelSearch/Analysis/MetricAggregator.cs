#nullable enable
using System;
using System.Collections.Generic;

namespace TaskLayer.ParallelSearch.Analysis;


/// <summary>
/// Aggregates results from multiple collectors into a single result object
/// compatible with ITransientDbResults for caching
/// </summary>
public class MetricAggregator
{
    private readonly List<IMetricCollector> _collectors = [];

    public MetricAggregator(IEnumerable<IMetricCollector> collectors)
    {
        _collectors.AddRange(collectors);
    }

    /// <summary>
    /// Runs all collectors on the context and produces an aggregated result
    /// </summary>
    public TransientDatabaseMetrics RunAnalysis(TransientDatabaseContext context)
    {
        var result = new TransientDatabaseMetrics(context.DatabaseName);

        foreach (var collector in _collectors)
        {
            try
            {
                if (!collector.CanCollectData(context))
                {
                    // Skip this collector or log warning
                    result.Errors.Add($"Skipping analyzer {collector.CollectorName} due to insufficient data.");
                    continue;
                }

                var analysisResults = collector.CollectData(context);

                // Merge results into the aggregated result
                foreach (var kvp in analysisResults)
                {
                    result.Results[kvp.Key] = kvp.Value /*is double.NaN ? 0 : kvp.Value*/;
                }
            }
            catch (Exception ex)
            {
                // Log error but continue with other collectors
                result.Errors.Add($"{collector.CollectorName}: {ex.Message}");
            }
        }

        // After running all collectors, populate the typed properties for CSV serialization
        result.PopulatePropertiesFromResults();

        return result;
    }
}
