#nullable enable
using System.Collections.Generic;

namespace TaskLayer.ParallelSearch.Analysis;

/// <summary>
/// Interface for post-hoc analysis of transient database search results
/// Each analyzer produces one or more named metrics that can be aggregated across databases
/// </summary>
public interface IMetricCollector
{
    /// <summary>
    /// Unique name for this analyzer (e.g., "OrganismSpecificity", "FdrMetrics")
    /// </summary>
    string CollectorName { get; }

    /// <summary>
    /// Returns the column headers this analyzer will produce
    /// </summary>
    IEnumerable<string> GetOutputColumns();

    /// <summary>
    /// Performs analysis on the context and returns key-value pairs for each metric
    /// Keys should match the output columns returned by GetOutputColumns()
    /// </summary>
    Dictionary<string, object> CollectData(TransientDatabaseContext context);

    /// <summary>
    /// Validates that all required data is present in the context
    /// Returns true if the analyzer can run, false otherwise
    /// </summary>
    bool CanCollectData(TransientDatabaseContext context);
}