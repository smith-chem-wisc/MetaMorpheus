#nullable enable
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
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

    // One writer formats every checkpoint row. CsvHelper compiles a write delegate for TransientDatabaseMetrics (166
    // columns) per writer, so a writer per row cost about 58 ms each, inside _writeLock, which held every parallel
    // search in line behind it. Both are used only under _writeLock.
    private static readonly UTF8Encoding Utf8NoBom = new(false);
    private StringWriter? _rowBuffer;
    private CsvWriter? _rowWriter;
    private string? _headerLine;

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
    /// A user-facing description of rows <see cref="InitializeCache"/> left out because older collectors wrote them
    /// (<see cref="TransientDatabaseMetrics.CurrentCollectorVersion"/>), or null when there were none.
    /// </summary>
    public string? OlderCollectorsMessage { get; private set; }

    /// <summary>Where the summary holding those rows was copied before they are dropped from it, or null.</summary>
    public string? OlderCollectorsBackupPath { get; private set; }

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

                // rows from older collectors are searched again, never reused (see CurrentCollectorVersion)
                var older = results.Where(r => r.CollectorVersion < TransientDatabaseMetrics.CurrentCollectorVersion).ToList();
                if (older.Count > 0)
                {
                    OlderCollectorsBackupPath = CopyAside("olderCollectors");
                    OlderCollectorsMessage = $"{older.Count} of {results.Count} rows in {Path.GetFileName(_csvFilePath)} were " +
                        $"written by older collectors (version {older.Min(r => r.CollectorVersion)}, current " +
                        $"{TransientDatabaseMetrics.CurrentCollectorVersion}); those databases will be searched again" +
                        (OlderCollectorsBackupPath is null ? "." : $". The original summary was copied to {OlderCollectorsBackupPath}.");
                    results = results.Except(older).ToList();
                }

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
    /// Copies the summary to the first free <c>name_label.csv</c> (then <c>name_label_N.csv</c>) beside it and returns the
    /// copy's path, or null when the copy failed.
    /// </summary>
    private string? CopyAside(string label)
    {
        try
        {
            string directory = Path.GetDirectoryName(_csvFilePath) ?? string.Empty;
            string baseName = Path.GetFileNameWithoutExtension(_csvFilePath);
            string extension = Path.GetExtension(_csvFilePath);
            string copyPath = Path.Combine(directory, $"{baseName}_{label}{extension}");
            for (int counter = 2; File.Exists(copyPath); counter++)
                copyPath = Path.Combine(directory, $"{baseName}_{label}_{counter}{extension}");
            File.Copy(_csvFilePath, copyPath);
            return copyPath;
        }
        catch (Exception)
        {
            return null;
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

            // formatted before the file is touched, so a row that cannot be formatted leaves the file as it was
            string text = rollbackPosition == 0 ? GetHeaderLine() + FormatRow(result) : FormatRow(result);
            byte[] bytes = Utf8NoBom.GetBytes(text);

            try
            {
                using (var stream = new FileStream(_csvFilePath, FileMode.Append, FileAccess.Write, FileShare.Read))
                {
                    stream.Write(bytes, 0, bytes.Length);
                }

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

    private static CsvConfiguration CreateWriteConfiguration(bool hasHeaderRecord) =>
        new(CultureInfo.InvariantCulture) { HasHeaderRecord = hasHeaderRecord };

    /// <summary>The header row, line ending included; built once, since a file that is deleted mid-run needs it again.</summary>
    private string GetHeaderLine()
    {
        if (_headerLine is null)
        {
            using var buffer = new StringWriter(CultureInfo.InvariantCulture);
            using (var csv = new CsvWriter(buffer, CreateWriteConfiguration(true), leaveOpen: true))
            {
                csv.WriteHeader<TransientDatabaseMetrics>();
                csv.NextRecord();
                csv.Flush();
            }
            _headerLine = buffer.ToString();
        }
        return _headerLine;
    }

    /// <summary>One data row, line ending included, written by the reused row writer.</summary>
    private string FormatRow(TransientDatabaseMetrics result)
    {
        if (_rowWriter is null)
        {
            _rowBuffer = new StringWriter(CultureInfo.InvariantCulture);
            _rowWriter = new CsvWriter(_rowBuffer, CreateWriteConfiguration(false), leaveOpen: true);
        }

        try
        {
            _rowWriter.WriteRecord(result);
            _rowWriter.NextRecord();
            _rowWriter.Flush();
            return _rowBuffer!.ToString();
        }
        catch
        {
            // a row that failed part way may have left fields in the writer's buffer; start the next row on a new writer
            _rowWriter.Dispose();
            _rowWriter = null;
            throw;
        }
        finally
        {
            _rowBuffer?.GetStringBuilder().Clear();
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
