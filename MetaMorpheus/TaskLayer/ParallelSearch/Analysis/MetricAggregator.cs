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
                // a collector without the data it needs (no protein groups without parsimony, say) leaves its columns
                // at their defaults; that is expected, not an error
                if (!collector.CanCollectData(context))
                {
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
                // keep going with the other collectors; the task warns about the failure and does not cache the row
                result.Errors.Add($"{collector.CollectorName}: {ex.Message}");
            }
        }

        // After running all collectors, populate the typed properties for CSV serialization
        result.PopulatePropertiesFromResults();

        return result;
    }
}
