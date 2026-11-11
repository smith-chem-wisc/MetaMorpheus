using System.IO;
using NUnit.Framework;
using TaskLayer.ParallelSearch;
using TaskLayer.ParallelSearch.Analysis;

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

    private static string GetTempCsvPath()
    {
        return Path.Combine(Path.GetTempPath(), $"parallel_search_cache_{Path.GetRandomFileName()}.csv");
    }
}
