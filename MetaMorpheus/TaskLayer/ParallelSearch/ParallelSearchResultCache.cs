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
    /// Where an existing summary that could not be read was moved by <see cref="InitializeCache"/>,
    /// or null when there was nothing to move or the move failed.
    /// </summary>
    public string? UnreadableFileBackupPath { get; private set; }

    /// <summary>
    /// A user-facing description of an existing summary that could not be read, or null when it was read
    /// (or did not exist).
    /// </summary>
    public string? UnreadableFileMessage { get; private set; }

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
        catch (Exception readException)
        {
            // Start fresh, but keep the unreadable file so the final rewrite does not destroy earlier results
            lock (_cacheLock)
            {
                _completedDatabases.Clear();
                _databaseResults.Clear();
            }

            BackUpUnreadableFile(readException);
        }
    }

    /// <summary>
    /// Moves the unreadable summary to the first free <c>name_N.csv</c> beside it (N = 1, 2, ...), so
    /// checkpoints start a fresh file and the old rows survive for the user to recover.
    /// </summary>
    private void BackUpUnreadableFile(Exception readException)
    {
        string fileName = Path.GetFileName(_csvFilePath);
        try
        {
            string directory = Path.GetDirectoryName(_csvFilePath) ?? string.Empty;
            string baseName = Path.GetFileNameWithoutExtension(_csvFilePath);
            string extension = Path.GetExtension(_csvFilePath);

            int counter = 1;
            string backupPath;
            do
            {
                backupPath = Path.Combine(directory, $"{baseName}_{counter}{extension}");
                counter++;
            } while (File.Exists(backupPath));

            lock (_writeLock)
            {
                File.Move(_csvFilePath, backupPath);
            }

            UnreadableFileBackupPath = backupPath;
            UnreadableFileMessage = $"Could not read the existing {fileName} ({readException.Message}). "
                + $"It was moved to {Path.GetFileName(backupPath)} and every transient database will be treated as not yet searched.";
        }
        catch (Exception backupException)
        {
            UnreadableFileBackupPath = null;
            UnreadableFileMessage = $"Could not read the existing {fileName} ({readException.Message}), "
                + $"and could not back it up ({backupException.Message}). It will be overwritten when the task finishes.";
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
