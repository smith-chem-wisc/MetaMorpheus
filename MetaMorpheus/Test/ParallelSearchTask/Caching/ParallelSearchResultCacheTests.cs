using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CsvHelper;
using CsvHelper.Configuration;
using NUnit.Framework;
using TaskLayer.ParallelSearch;
using TaskLayer.ParallelSearch.Analysis;
using TaskLayer.ParallelSearch.Analysis.Collectors;

namespace Test.ParallelSearchTask.Caching;

[TestFixture]
public class ParallelSearchResultCacheTests
{
    [Test]
    public void Add_HoldsResultInMemoryWithoutMarkingItCompleted()
    {
        string path = GetTempCsvPath();
        var cache = new ParallelSearchResultCache(path);
        var metrics = new TransientDatabaseMetrics("Db1");

        bool firstAdd = cache.Add(metrics);
        bool secondAdd = cache.Add(metrics);
        bool found = cache.TryGetValue("Db1", out var held);
        bool removed = cache.Remove(metrics);

        // A database only counts as completed once its row is checkpointed to the file
        Assert.Multiple(() =>
        {
            Assert.That(firstAdd, Is.True);
            Assert.That(secondAdd, Is.True);
            Assert.That(found, Is.True);
            Assert.That(held, Is.SameAs(metrics));
            Assert.That(cache.Contains("Db1"), Is.False);
            Assert.That(cache.Count, Is.EqualTo(0));
            Assert.That(removed, Is.False);
            Assert.That(File.Exists(path), Is.False);
        });
    }

    [Test]
    public void AddAndWrite_ThenRemove_TracksCompletedDatabases()
    {
        string path = GetTempCsvPath();
        try
        {
            var cache = new ParallelSearchResultCache(path);
            var metrics = new TransientDatabaseMetrics("Db1");

            bool written = cache.AddAndWrite(metrics);
            int countAfterWrite = cache.Count;
            bool containsAfterWrite = cache.Contains("Db1");
            bool removed = cache.Remove(metrics);

            Assert.Multiple(() =>
            {
                Assert.That(written, Is.True);
                Assert.That(countAfterWrite, Is.EqualTo(1));
                Assert.That(containsAfterWrite, Is.True);
                Assert.That(removed, Is.True);
                Assert.That(cache.Count, Is.EqualTo(0));
                Assert.That(cache.TryGetValue("Db1", out _), Is.False);
            });
        }
        finally
        {
            if (File.Exists(path))
            {
                File.Delete(path);
            }
        }
    }

    [Test]
    public void InitializeCache_LeavesOutRowsFromOlderCollectors_AndCopiesTheSummaryAside()
    {
        string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, $"parallel_search_cache_{Path.GetRandomFileName()}");
        Directory.CreateDirectory(folder);
        string path = Path.Combine(folder, TransientDatabaseResultsManager.SummaryResultsFileName);
        try
        {
            var writer = new ParallelSearchResultCache(path);
            Assert.That(writer.AddAndWrite(new TransientDatabaseMetrics("Current")), Is.True);
            Assert.That(writer.AddAndWrite(new TransientDatabaseMetrics("Older") { CollectorVersion = 0 }), Is.True);

            var cache = new ParallelSearchResultCache(path);
            cache.InitializeCache();

            Assert.Multiple(() =>
            {
                Assert.That(new TransientDatabaseMetrics("New").CollectorVersion, Is.EqualTo(TransientDatabaseMetrics.CurrentCollectorVersion));
                Assert.That(cache.Contains("Current"), Is.True);
                Assert.That(cache.Contains("Older"), Is.False, "a row from older collectors is searched again");
                Assert.That(cache.TryGetValue("Current", out var loaded) ? loaded!.CollectorVersion : -1,
                    Is.EqualTo(TransientDatabaseMetrics.CurrentCollectorVersion));
                Assert.That(cache.OlderCollectorsMessage, Does.Contain("1 of 2 rows"));
                Assert.That(cache.OlderCollectorsBackupPath, Is.EqualTo(Path.Combine(folder, "ManySearchSummary_olderCollectors.csv")));
                Assert.That(File.Exists(cache.OlderCollectorsBackupPath), Is.True);
            });
        }
        finally
        {
            Directory.Delete(folder, true);
        }
    }

    [Test]
    public void InitializeCache_ReadsASummaryWithoutTheVersionColumnAsOlder()
    {
        string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, $"parallel_search_cache_{Path.GetRandomFileName()}");
        Directory.CreateDirectory(folder);
        string path = Path.Combine(folder, TransientDatabaseResultsManager.SummaryResultsFileName);
        try
        {
            // a summary written before the column existed: write one and drop the column
            var writer = new ParallelSearchResultCache(path);
            Assert.That(writer.AddAndWrite(new TransientDatabaseMetrics("OldDb")), Is.True);
            string[] lines = File.ReadAllLines(path);
            int column = System.Array.IndexOf(lines[0].Split(','), nameof(TransientDatabaseMetrics.CollectorVersion));
            Assert.That(column, Is.GreaterThanOrEqualTo(0));
            File.WriteAllLines(path, System.Linq.Enumerable.Select(lines,
                l => string.Join(",", System.Linq.Enumerable.Where(l.Split(','), (_, i) => i != column))));

            var cache = new ParallelSearchResultCache(path);
            cache.InitializeCache();

            Assert.Multiple(() =>
            {
                Assert.That(cache.Contains("OldDb"), Is.False);
                Assert.That(cache.OlderCollectorsMessage, Does.Contain("version 0"));
                Assert.That(cache.UnreadableFileMessage, Is.Null);
            });
        }
        finally
        {
            Directory.Delete(folder, true);
        }
    }

    [Test]
    public void AddAndWrite_ThenInitializeCache_LoadsPersistedResult()
    {
        string path = GetTempCsvPath();
        try
        {
            var firstCache = new ParallelSearchResultCache(path);
            var metrics = new TransientDatabaseMetrics("DbPersisted");

            Assert.That(firstCache.AddAndWrite(metrics), Is.True);

            var secondCache = new ParallelSearchResultCache(path);
            secondCache.InitializeCache();

            bool found = secondCache.TryGetValue("DbPersisted", out var loaded);

            Assert.Multiple(() =>
            {
                Assert.That(found, Is.True);
                Assert.That(loaded, Is.Not.Null);
                Assert.That(loaded!.DatabaseName, Is.EqualTo("DbPersisted"));
                Assert.That(secondCache.Count, Is.EqualTo(1));
            });
        }
        finally
        {
            if (File.Exists(path))
            {
                File.Delete(path);
            }
        }
    }

    [Test]
    public void InitializeCache_WhenFileCannotBeRead_ClearsInMemoryCache()
    {
        string path = GetTempCsvPath();
        try
        {
            var cache = new ParallelSearchResultCache(path);
            cache.AddAndWrite(new TransientDatabaseMetrics("DbInMemory"));
            Assert.That(cache.Count, Is.EqualTo(1));

            using (new FileStream(path, FileMode.Open, FileAccess.ReadWrite, FileShare.None))
            {
                cache.InitializeCache();
            }

            Assert.Multiple(() =>
            {
                Assert.That(cache.Count, Is.EqualTo(0));
                Assert.That(cache.TryGetValue("DbInMemory", out _), Is.False);
                Assert.That(cache.UnreadableFileBackupPath, Is.Null);
                Assert.That(cache.UnreadableFileMessage, Does.Contain("could not back it up"));
            });
        }
        finally
        {
            if (File.Exists(path))
            {
                File.Delete(path);
            }
        }
    }

    [Test]
    public void InitializeCache_WhenFileCannotBeParsed_MovesItToTheFirstFreeNumberedBackup()
    {
        string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, $"parallel_search_cache_{Path.GetRandomFileName()}");
        Directory.CreateDirectory(folder);
        string path = Path.Combine(folder, TransientDatabaseResultsManager.SummaryResultsFileName);
        try
        {
            byte[] firstBytes = Encoding.UTF8.GetBytes("this is not,a summary\nfile,at all\n");
            File.WriteAllBytes(path, firstBytes);

            var firstCache = new ParallelSearchResultCache(path);
            firstCache.InitializeCache();

            string firstBackup = Path.Combine(folder, "ManySearchSummary_1.csv");
            Assert.Multiple(() =>
            {
                Assert.That(firstCache.Count, Is.EqualTo(0));
                Assert.That(firstCache.UnreadableFileBackupPath, Is.EqualTo(firstBackup));
                Assert.That(firstCache.UnreadableFileMessage, Does.Contain("ManySearchSummary_1.csv"));
                Assert.That(File.ReadAllBytes(firstBackup), Is.EqualTo(firstBytes));
                Assert.That(File.Exists(path), Is.False);
            });

            byte[] secondBytes = Encoding.UTF8.GetBytes("another,broken\nsummary\n");
            File.WriteAllBytes(path, secondBytes);

            var secondCache = new ParallelSearchResultCache(path);
            secondCache.InitializeCache();

            string secondBackup = Path.Combine(folder, "ManySearchSummary_2.csv");
            Assert.Multiple(() =>
            {
                Assert.That(secondCache.UnreadableFileBackupPath, Is.EqualTo(secondBackup));
                Assert.That(File.ReadAllBytes(secondBackup), Is.EqualTo(secondBytes));
                Assert.That(File.ReadAllBytes(firstBackup), Is.EqualTo(firstBytes));
            });

            // a fresh run checkpoints into a new, readable file
            secondCache.AddAndWrite(new TransientDatabaseMetrics("DbAfterBackup"));
            var thirdCache = new ParallelSearchResultCache(path);
            thirdCache.InitializeCache();
            Assert.Multiple(() =>
            {
                Assert.That(thirdCache.UnreadableFileMessage, Is.Null);
                Assert.That(thirdCache.Contains("DbAfterBackup"), Is.True);
            });
        }
        finally
        {
            if (Directory.Exists(folder))
            {
                Directory.Delete(folder, true);
            }
        }
    }

    [Test]
    public void InitializeCache_WhenFileIsReadable_ReportsNothing()
    {
        string path = GetTempCsvPath();
        try
        {
            new ParallelSearchResultCache(path).AddAndWrite(new TransientDatabaseMetrics("DbReadable"));

            var cache = new ParallelSearchResultCache(path);
            cache.InitializeCache();

            Assert.Multiple(() =>
            {
                Assert.That(cache.UnreadableFileMessage, Is.Null);
                Assert.That(cache.UnreadableFileBackupPath, Is.Null);
                Assert.That(cache.Contains("DbReadable"), Is.True);
            });
        }
        finally
        {
            if (File.Exists(path))
            {
                File.Delete(path);
            }
        }
    }

    [Test]
    public void AddAndWrite_ThenInitializeCache_KeepsTheRetentionTimeCorrelationCoefficients()
    {
        string path = GetTempCsvPath();
        try
        {
            var metrics = new TransientDatabaseMetrics("DbRt");
            metrics.Results[RetentionTimeCollector.PsmRtCorrelationCoefficient] = 0.91;
            metrics.Results[RetentionTimeCollector.PeptideRtCorrelationCoefficient] = 0.87;
            Assert.That(new ParallelSearchResultCache(path).AddAndWrite(metrics), Is.True);

            var cache = new ParallelSearchResultCache(path);
            cache.InitializeCache();
            Assert.That(cache.TryGetValue("DbRt", out var loaded), Is.True);

            Assert.Multiple(() =>
            {
                Assert.That(loaded!.Psm_RtCorrelationCoefficient, Is.EqualTo(0.91).Within(1e-9));
                Assert.That(loaded.Peptide_RtCorrelationCoefficient, Is.EqualTo(0.87).Within(1e-9));
                Assert.That((double)loaded.Results[RetentionTimeCollector.PsmRtCorrelationCoefficient], Is.EqualTo(0.91).Within(1e-9));
                Assert.That((double)loaded.Results[RetentionTimeCollector.PeptideRtCorrelationCoefficient], Is.EqualTo(0.87).Within(1e-9));
            });
        }
        finally
        {
            if (File.Exists(path))
            {
                File.Delete(path);
            }
        }
    }

    [Test]
    public void AppendToFile_WritesTheSameBytesAsAWriterPerRow()
    {
        string path = GetTempCsvPath();
        string expectedPath = GetTempCsvPath();
        try
        {
            var cache = new ParallelSearchResultCache(path);
            List<TransientDatabaseMetrics> rows = Enumerable.Range(0, 5).Select(MakeMetrics).ToList();
            foreach (var row in rows)
            {
                cache.AppendToFile(row);
            }

            // the writer-per-row checkpoint this replaced, kept here as the reference for the file's bytes
            foreach (var row in rows)
            {
                AppendWithAWriterPerRow(expectedPath, row);
            }

            Assert.That(File.ReadAllBytes(path), Is.EqualTo(File.ReadAllBytes(expectedPath)));
        }
        finally
        {
            DeleteIfPresent(path);
            DeleteIfPresent(expectedPath);
        }
    }

    [Test]
    public void AppendToFile_ManyRows_WritesOneHeaderAndReadsBackEveryRow()
    {
        string path = GetTempCsvPath();
        try
        {
            var cache = new ParallelSearchResultCache(path);
            for (int i = 0; i < 50; i++)
            {
                cache.AppendToFile(MakeMetrics(i));
            }

            string header = File.ReadLines(path).First();
            var reread = new ParallelSearchResultCache(path);
            reread.InitializeCache();

            Assert.Multiple(() =>
            {
                Assert.That(File.ReadLines(path).Count(line => line == header), Is.EqualTo(1));
                Assert.That(reread.Count, Is.EqualTo(50));
                Assert.That(reread.TryGetValue("Db, \"49\"", out var last), Is.True);
                Assert.That(last!.Psm_RtCorrelationCoefficient, Is.EqualTo(0.5 + 49 / 1000.0).Within(1e-9));
            });
        }
        finally
        {
            DeleteIfPresent(path);
        }
    }

    [Test]
    public void AppendToFile_AfterTheFileIsDeleted_StartsTheNewFileWithTheHeader()
    {
        string path = GetTempCsvPath();
        try
        {
            var cache = new ParallelSearchResultCache(path);
            cache.AppendToFile(MakeMetrics(0));
            string header = File.ReadLines(path).First();

            File.Delete(path);
            cache.AppendToFile(MakeMetrics(1));

            var reread = new ParallelSearchResultCache(path);
            reread.InitializeCache();
            Assert.Multiple(() =>
            {
                Assert.That(File.ReadLines(path).First(), Is.EqualTo(header));
                Assert.That(reread.Count, Is.EqualTo(1));
                Assert.That(reread.Contains("Db, \"1\""), Is.True);
            });
        }
        finally
        {
            DeleteIfPresent(path);
        }
    }

    [Test]
    public void AppendToFile_WhenTheFileCannotBeOpened_ThrowsAndLeavesItAsItWas()
    {
        string path = GetTempCsvPath();
        try
        {
            var cache = new ParallelSearchResultCache(path);
            cache.AppendToFile(MakeMetrics(0));
            byte[] before = File.ReadAllBytes(path);

            using (new FileStream(path, FileMode.Open, FileAccess.ReadWrite, FileShare.None))
            {
                Assert.That(() => cache.AppendToFile(MakeMetrics(1)), Throws.InstanceOf<IOException>());
            }

            // the failed row is not counted, and the next row still writes cleanly on the same cache
            Assert.That(cache.Contains("Db, \"1\""), Is.False);
            Assert.That(File.ReadAllBytes(path), Is.EqualTo(before));
            cache.AppendToFile(MakeMetrics(2));

            var reread = new ParallelSearchResultCache(path);
            reread.InitializeCache();
            Assert.Multiple(() =>
            {
                Assert.That(reread.Count, Is.EqualTo(2));
                Assert.That(reread.Contains("Db, \"2\""), Is.True);
            });
        }
        finally
        {
            DeleteIfPresent(path);
        }
    }

    [Test]
    public void AppendToFile_FromManyThreads_WritesEveryRowWhole()
    {
        string path = GetTempCsvPath();
        try
        {
            var cache = new ParallelSearchResultCache(path);
            Parallel.For(0, 400, new ParallelOptions { MaxDegreeOfParallelism = 16 }, i => cache.AppendToFile(MakeMetrics(i)));

            var reread = new ParallelSearchResultCache(path);
            reread.InitializeCache();
            Assert.Multiple(() =>
            {
                Assert.That(cache.Count, Is.EqualTo(400));
                Assert.That(reread.Count, Is.EqualTo(400));
                for (int i = 0; i < 400; i++)
                {
                    Assert.That(reread.TryGetValue($"Db, \"{i}\"", out var row), Is.True, $"row {i}");
                    Assert.That(row!.Psm_RtCorrelationCoefficient, Is.EqualTo(0.5 + i / 1000.0).Within(1e-9), $"row {i}");
                }
            });
        }
        finally
        {
            DeleteIfPresent(path);
        }
    }

    /// <summary>A row whose name needs quoting and whose values differ per row.</summary>
    private static TransientDatabaseMetrics MakeMetrics(int i)
    {
        var metrics = new TransientDatabaseMetrics($"Db, \"{i}\"");
        metrics.Results[RetentionTimeCollector.PsmRtCorrelationCoefficient] = 0.5 + i / 1000.0;
        metrics.Results[RetentionTimeCollector.PeptideRtCorrelationCoefficient] = 0.25 + i / 1000.0;
        metrics.AllPsmsPerProteinGroup = new[] { i, i + 0.5, 3.0 };
        return metrics;
    }

    private static void AppendWithAWriterPerRow(string path, TransientDatabaseMetrics result)
    {
        result.PopulatePropertiesFromResults();
        bool isNew = !File.Exists(path) || new FileInfo(path).Length == 0;
        using var stream = new FileStream(path, FileMode.Append, FileAccess.Write, FileShare.Read);
        using var writer = new StreamWriter(stream);
        using var csv = new CsvWriter(writer, new CsvConfiguration(CultureInfo.InvariantCulture) { HasHeaderRecord = isNew });
        if (isNew)
        {
            csv.WriteHeader<TransientDatabaseMetrics>();
            csv.NextRecord();
        }
        csv.WriteRecord(result);
        csv.NextRecord();
        csv.Flush();
    }

    private static void DeleteIfPresent(string path)
    {
        if (File.Exists(path))
        {
            File.Delete(path);
        }
    }

    private static string GetTempCsvPath()
    {
        return Path.Combine(TestContext.CurrentContext.TestDirectory, $"parallel_search_cache_{Path.GetRandomFileName()}.csv");
    }
}
