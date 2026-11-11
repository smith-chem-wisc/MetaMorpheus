#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using TaskLayer.ParallelSearch.Analysis;

namespace TaskLayer.ParallelSearch;

/// <summary>
/// Runs the metric collectors on each transient database and keeps the results in the summary file,
/// so a rerun can skip databases that already have a row.
/// </summary>
public class TransientDatabaseResultsManager
{
    private readonly MetricAggregator _metricAggregator;
    private readonly ParallelSearchResultCache _analysisCache;

    public const string SummaryResultsFileName = "ManySearchSummary.csv";

    /// <summary>
    /// Gets the number of databases with cached analysis results
    /// </summary>
    public int CachedAnalysisCount => _analysisCache.Count;
    public string SummaryFilePath => _analysisCache.FilePath;
    
    /// <summary>
    /// Gets all cached analysis results
    /// </summary>
    public Dictionary<string, TransientDatabaseMetrics> TransientDatabaseMetricsDictionary => _analysisCache.AllResultsDictionary;
    public void PersistAnalysisCache()
    {
        _analysisCache.WriteAllToFile();
    }

    public void AppendCheckpoint(TransientDatabaseMetrics result)
    {
        _analysisCache.AppendToFile(result);
    }

    /// <param name="metricAggregator">Aggregator for computing metrics (counts, organism specificity, etc.)</param>
    /// <param name="analysisCachePath">Path to the summary CSV, which doubles as the cache of finished databases</param>
    public TransientDatabaseResultsManager(MetricAggregator metricAggregator, string analysisCachePath)
    {
        _metricAggregator = metricAggregator ?? throw new ArgumentNullException(nameof(metricAggregator));
        _analysisCache = new ParallelSearchResultCache(analysisCachePath);
        _analysisCache.InitializeCache();
    }

    #region Search Result Analysis Cache Methods

    /// <summary>
    /// Checks if a database has cached analysis results
    /// </summary>
    /// <param name="databaseName">Name of the database to check</param>
    /// <returns>True if results are cached, false otherwise</returns>
    public bool HasCachedResults(string databaseName)
    {
        return _analysisCache.Contains(databaseName);
    }

    /// <summary>
    /// Attempts to retrieve cached analysis results for a database
    /// </summary>
    /// <param name="databaseName">Name of the database</param>
    /// <param name="result">Cached result if found, null otherwise</param>
    /// <returns>True if results were found in cache</returns>
    public bool TryGetCachedResult(string databaseName, out TransientDatabaseMetrics? result)
    {
        return _analysisCache.TryGetValue(databaseName, out result);
    }

    /// <summary>
    /// Removes a database's results from cache
    /// Use this when overwriting results for a specific database
    /// </summary>
    /// <param name="databaseName">Name of the database to remove</param>
    /// <returns>True if the database was in cache and removed</returns>
    public bool RemoveFromCache(string databaseName)
    {
        if (!_analysisCache.TryGetValue(databaseName, out var result) || result == null)
            return false;

        return _analysisCache.Remove(result);
    }

    /// <summary>
    /// Processes a single database: runs analysis and caches results
    /// This method is thread-safe and can be called in parallel for multiple databases
    /// </summary>
    /// <param name="context">Analysis context containing PSMs, peptides, and protein groups</param>
    /// <param name="forceRecompute">If true, recomputes even if cached (removes old cache entry)</param>
    /// <returns>Analysis results (either cached or newly computed)</returns>
    public TransientDatabaseMetrics ProcessDatabase(
        TransientDatabaseContext context,
        bool forceRecompute = false)
    {
        if (context == null)
            throw new ArgumentNullException(nameof(context));

        string dbName = context.DatabaseName;

        // Check cache first (unless forcing recompute)
        if (!forceRecompute && TryGetCachedResult(dbName, out var cachedResult) && cachedResult != null)
        {
            return cachedResult;
        }

        // If forcing recompute, remove old cache entry
        if (forceRecompute)
        {
            RemoveFromCache(dbName);
        }

        // Run analysis
        var analysisResult = _metricAggregator.RunAnalysis(context);

        // Cache analysis results in memory immediately; durable checkpointing happens
        // after per-database output files are written successfully.
        _analysisCache.Add(analysisResult);

        return analysisResult;
    }

    /// <summary>
    /// Gets summary statistics about cache status
    /// </summary>
    public CacheSummary GetCacheSummary(List<string> allDatabaseNames)
    {
        int cached = allDatabaseNames.Count(HasCachedResults);
        int needsProcessing = allDatabaseNames.Count - cached;

        return new CacheSummary
        {
            TotalDatabases = allDatabaseNames.Count,
            CachedDatabases = cached,
            DatabasesNeedingProcessing = needsProcessing,
            CacheHitRate = allDatabaseNames.Count > 0
                ? (double)cached / allDatabaseNames.Count
                : 0.0
        };
    }

    #endregion
}

/// <summary>
/// Summary of cache status for reporting
/// </summary>
public class CacheSummary
{
    public int TotalDatabases { get; init; }
    public int CachedDatabases { get; init; }
    public int DatabasesNeedingProcessing { get; init; }
    public double CacheHitRate { get; init; }

    public override string ToString()
    {
        return $"Cache Status: {CachedDatabases}/{TotalDatabases} databases cached " +
               $"({CacheHitRate:P1} hit rate), {DatabasesNeedingProcessing} need processing";
    }
}
