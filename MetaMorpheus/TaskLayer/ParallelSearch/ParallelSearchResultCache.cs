#nullable enable
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using CsvHelper;
using CsvHelper.Configuration;
using TaskLayer.ParallelSearch.Analysis;

namespace TaskLayer.ParallelSearch;

/// <summary>
/// Helper class for thread-safe reading and writing of database search results to CSV
/// Manages caching, validation, and tracking of database search results
/// </summary>
public class ParallelSearchResultCache
{
    private readonly object _writeLock = new();
    private readonly object _cacheLock = new();
    private readonly HashSet<string> _completedDatabases = new();
    private readonly ConcurrentDictionary<string, TransientDatabaseMetrics> _databaseResults = new();
    private readonly string _csvFilePath;

    public string FilePath => _csvFilePath;

    /// <summary>
    /// Helper class for thread-safe reading and writing of database search results to CSV
    /// Manages caching, validation, and tracking of database search results
    /// </summary>
    public ParallelSearchResultCache(string csvFilePath)
    {
        _csvFilePath = csvFilePath;
    }

    /// <summary>
    /// Gets the number of completed databases
    /// </summary>
    public int Count
    {
        get
        {
            lock (_cacheLock)
            {
                return _completedDatabases.Count;
            }
        }
    }

    /// <summary>
    /// Gets all results currently in cache
    /// </summary>
    public Dictionary<string, TransientDatabaseMetrics> AllResultsDictionary => _databaseResults.ToDictionary(kvp => kvp.Key, kvp => kvp.Value);

    public List<TransientDatabaseMetrics> AllResultsList => _databaseResults.Values.ToList();


    #region Dictionary Like Methods

    public bool Add(TransientDatabaseMetrics? result)
    {
        if (result is null) return false;

        string key = result.DatabaseName;


        lock (_cacheLock)
        {
            _databaseResults[key] = result;
        }
        return true;
    }

    public bool Remove(TransientDatabaseMetrics? result)
    {
        if (result is null) return false;

        string key = result.DatabaseName;
        if (!_databaseResults.ContainsKey(key))
            return false;

        lock (_cacheLock)
        {
            if (!_completedDatabases.Contains(key))
                return false;

            _completedDatabases.Remove(key);
            _databaseResults.Remove(key, out _);
        }
        return true;
    }

    public bool AddAndWrite(TransientDatabaseMetrics? result)
    {
        if (result is null) return false;

        if (!Add(result))
            return false;

        AppendToFile(result);
        return true;
    }

    public bool TryGetValue(string key, out TransientDatabaseMetrics? result) => _databaseResults.TryGetValue(key, out result);

    /// <summary>
    /// Checks if a database result already exists in cache
    /// </summary>
    public bool Contains(string databaseName)
    {
        lock (_cacheLock)
        {
            return _completedDatabases.Contains(databaseName);
        }
    }

    #endregion

    #region IO 

    /// <summary>
    /// Loads cached results from the CSV file if it exists and returns both the list and count
    /// </summary>
    public void InitializeCache()
    {
        if (!File.Exists(_csvFilePath))
            return;

        try
        {
            lock (_writeLock)
            {
                using var reader = new StreamReader(_csvFilePath);
                using var csv = new CsvReader(reader, new CsvConfiguration(CultureInfo.InvariantCulture)
                {
                    HasHeaderRecord = true,
                    MissingFieldFound = null,
                    ReadingExceptionOccurred = args => false
                });

                var results = csv.GetRecords<TransientDatabaseMetrics>().ToList();

                // Populate Results dictionary from properties
                foreach (var result in results)
                {
                    result.PopulateResultsFromProperties();
                }

                lock (_cacheLock)
                {
                    foreach (var result in results)
                    {
                        _completedDatabases.Add(result.DatabaseName);
                        _databaseResults[result.DatabaseName] = result;
                    }
                }
            }
        }
        catch
        {
            // If there's an error reading the cache, start fresh
            lock (_cacheLock)
            {
                _completedDatabases.Clear();
                _databaseResults.Clear();
            }
        }
    }

    /// <summary>
    /// Writes a single result to the CSV file in a thread-safe manner and tracks it
    /// </summary>
    public void AppendToFile(TransientDatabaseMetrics result)
    {
        lock (_writeLock)
        {
            result.PopulatePropertiesFromResults();

            long rollbackPosition = 0;
            if (File.Exists(_csvFilePath))
            {
                rollbackPosition = new FileInfo(_csvFilePath).Length;
            }

            try
            {
                using var stream = new FileStream(_csvFilePath, FileMode.Append, FileAccess.Write, FileShare.Read);
                using var writer = new StreamWriter(stream);
                using var csv = new CsvWriter(writer, new CsvConfiguration(CultureInfo.InvariantCulture)
                {
                    HasHeaderRecord = rollbackPosition == 0
                });

                // Write header if file is new
                if (rollbackPosition == 0)
                {
                    csv.WriteHeader<TransientDatabaseMetrics>();
                    csv.NextRecord();
                }

                csv.WriteRecord(result);
                csv.NextRecord();
                csv.Flush();

                lock (_cacheLock)
                {
                    _completedDatabases.Add(result.DatabaseName);
                }
            }
            catch
            {
                // Roll back to pre-write position on failure to prevent partial rows
                if (rollbackPosition > 0)
                {
                    using var fs = new FileStream(_csvFilePath, FileMode.Open, FileAccess.Write);
                    fs.SetLength(rollbackPosition);
                }
                else if (File.Exists(_csvFilePath))
                {
                    File.Delete(_csvFilePath);
                }
                throw;
            }
        }
    }

    public void WriteAllToFile(string? outputPath = null)
    {
        outputPath ??= _csvFilePath;
        lock (_writeLock)
        {
            using var writer = new StreamWriter(outputPath);
            using var csv = new CsvWriter(writer, new CsvConfiguration(CultureInfo.InvariantCulture)
            {
                HasHeaderRecord = true
            });
            csv.WriteHeader<TransientDatabaseMetrics>();
            csv.NextRecord();
            foreach (var result in _databaseResults.Values
                .OrderByDescending(p => p.TargetPsmsFromTransientDbAtQValueThreshold)
                .ThenBy(p => p.DatabaseName, StringComparer.Ordinal))
            {
                result.PopulatePropertiesFromResults();

                csv.WriteRecord(result);
                csv.NextRecord();
            }
            csv.Flush();
        }
    }

    #endregion
}
