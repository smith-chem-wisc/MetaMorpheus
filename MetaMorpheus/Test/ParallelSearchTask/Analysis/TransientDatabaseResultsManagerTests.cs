using System;
using System.Collections.Generic;
using System.IO;
using NUnit.Framework;
using TaskLayer.ParallelSearch;

namespace Test.ParallelSearchTask.Analysis;

public class TransientDatabaseResultsManagerTests
{
    private string _testOutputDir = null!;

    [SetUp]
    public void SetUp()
    {
        _testOutputDir = Path.Combine(Path.GetTempPath(), $"ParallelSearchResults_{Guid.NewGuid():N}");
        Directory.CreateDirectory(_testOutputDir);
    }

    [TearDown]
    public void TearDown()
    {
        if (Directory.Exists(_testOutputDir))
            Directory.Delete(_testOutputDir, true);
    }

    [Test]
    public void Constructor_WithNullAggregator_Throws()
    {
        var ex = Assert.Throws<ArgumentNullException>(() => new TransientDatabaseResultsManager(
            null!,
            Path.Combine(_testOutputDir, "cache.csv")));

        Assert.That(ex.ParamName, Is.EqualTo("metricAggregator"));
    }

    [Test]
    public void AppendCheckpoint_LeavesOutARowWithCollectorErrors_SoARerunRetriesIt()
    {
        var aggregator = new TaskLayer.ParallelSearch.Analysis.MetricAggregator(new List<TaskLayer.ParallelSearch.Analysis.IMetricCollector>());
        string summaryPath = Path.Combine(_testOutputDir, "cache.csv");
        var manager = new TransientDatabaseResultsManager(aggregator, summaryPath);
        var failed = new TaskLayer.ParallelSearch.Analysis.TransientDatabaseMetrics("Failed");
        failed.Errors.Add("Thrower: Expected test failure");
        var succeeded = new TaskLayer.ParallelSearch.Analysis.TransientDatabaseMetrics("Succeeded");

        Assert.Multiple(() =>
        {
            Assert.That(manager.AppendCheckpoint(failed), Is.False);
            Assert.That(manager.AppendCheckpoint(succeeded), Is.True);
            Assert.That(manager.HasCachedResults("Failed"), Is.False);
            Assert.That(manager.HasCachedResults("Succeeded"), Is.True);
            Assert.That(new TransientDatabaseResultsManager(aggregator, summaryPath).HasCachedResults("Failed"), Is.False);
        });
    }

    [Test]
    public void ProcessDatabase_WhenACollectorFails_TheRowIsLeftOutOfTheSummaryWrittenAtTheEnd()
    {
        var aggregator = new TaskLayer.ParallelSearch.Analysis.MetricAggregator(
            new List<TaskLayer.ParallelSearch.Analysis.IMetricCollector> { new ThrowingCollector() });
        string summaryPath = Path.Combine(_testOutputDir, "cache.csv");
        var manager = new TransientDatabaseResultsManager(aggregator, summaryPath);
        var commonParameters = Utility.ParallelSearchTestContextFactory.CreateCommonParameters();
        var context = Utility.ParallelSearchTestContextFactory.CreateContext(commonParameters, allPsms: [], transientPsms: [],
            allPeptides: [], transientPeptides: []);

        var result = manager.ProcessDatabase(context);
        bool appended = manager.AppendCheckpoint(result);
        manager.PersistAnalysisCache();

        Assert.Multiple(() =>
        {
            Assert.That(result.Errors, Is.Not.Empty);
            Assert.That(appended, Is.False);
            Assert.That(new TransientDatabaseResultsManager(aggregator, summaryPath).HasCachedResults(context.DatabaseName), Is.False);
        });
    }

    [Test]
    public void HasCachedResults_WithNoCache_ReturnsFalse()
    {
        var aggregator = new TaskLayer.ParallelSearch.Analysis.MetricAggregator(new List<TaskLayer.ParallelSearch.Analysis.IMetricCollector>());
        var manager = new TransientDatabaseResultsManager(
            aggregator,
            Path.Combine(_testOutputDir, "cache.csv"));

        bool hasCache = manager.HasCachedResults("TestDb");

        Assert.That(hasCache, Is.False);
    }

    [Test]
    public void TryGetCachedResult_WithNoCache_ReturnsFalse()
    {
        var aggregator = new TaskLayer.ParallelSearch.Analysis.MetricAggregator(new List<TaskLayer.ParallelSearch.Analysis.IMetricCollector>());
        var manager = new TransientDatabaseResultsManager(
            aggregator,
            Path.Combine(_testOutputDir, "cache.csv"));

        bool found = manager.TryGetCachedResult("TestDb", out var result);

        Assert.Multiple(() =>
        {
            Assert.That(found, Is.False);
            Assert.That(result, Is.Null);
        });
    }

    [Test]
    public void ProcessDatabase_WithNullContext_Throws()
    {
        var aggregator = new TaskLayer.ParallelSearch.Analysis.MetricAggregator(new List<TaskLayer.ParallelSearch.Analysis.IMetricCollector>());
        var manager = new TransientDatabaseResultsManager(
            aggregator,
            Path.Combine(_testOutputDir, "cache.csv"));

        var ex = Assert.Throws<ArgumentNullException>(() => manager.ProcessDatabase(null!));

        Assert.That(ex.ParamName, Is.EqualTo("context"));
    }

    [Test]
    public void GetCacheSummary_EmptyList_ReturnsZeroCounts()
    {
        var aggregator = new TaskLayer.ParallelSearch.Analysis.MetricAggregator(new List<TaskLayer.ParallelSearch.Analysis.IMetricCollector>());
        var manager = new TransientDatabaseResultsManager(
            aggregator,
            Path.Combine(_testOutputDir, "cache.csv"));

        var summary = manager.GetCacheSummary(new List<string>());

        Assert.Multiple(() =>
        {
            Assert.That(summary.TotalDatabases, Is.EqualTo(0));
            Assert.That(summary.CachedDatabases, Is.EqualTo(0));
            Assert.That(summary.DatabasesNeedingProcessing, Is.EqualTo(0));
            Assert.That(summary.CacheHitRate, Is.EqualTo(0.0));
        });
    }

    [Test]
    public void GetCacheSummary_MultipleDatabases_ReturnsCorrectCounts()
    {
        var aggregator = new TaskLayer.ParallelSearch.Analysis.MetricAggregator(new List<TaskLayer.ParallelSearch.Analysis.IMetricCollector>());
        var manager = new TransientDatabaseResultsManager(
            aggregator,
            Path.Combine(_testOutputDir, "cache.csv"));

        var summary = manager.GetCacheSummary(new List<string> { "Db1", "Db2", "Db3" });

        Assert.Multiple(() =>
        {
            Assert.That(summary.TotalDatabases, Is.EqualTo(3));
            Assert.That(summary.CachedDatabases, Is.EqualTo(0));
            Assert.That(summary.DatabasesNeedingProcessing, Is.EqualTo(3));
            Assert.That(summary.CacheHitRate, Is.EqualTo(0.0));
        });
    }

    [Test]
    public void CachedAnalysisCount_InitiallyZero()
    {
        var aggregator = new TaskLayer.ParallelSearch.Analysis.MetricAggregator(new List<TaskLayer.ParallelSearch.Analysis.IMetricCollector>());
        var manager = new TransientDatabaseResultsManager(
            aggregator,
            Path.Combine(_testOutputDir, "cache.csv"));

        Assert.That(manager.CachedAnalysisCount, Is.EqualTo(0));
    }

    [Test]
    public void RemoveFromCache_WithNoCache_ReturnsFalse()
    {
        var aggregator = new TaskLayer.ParallelSearch.Analysis.MetricAggregator(new List<TaskLayer.ParallelSearch.Analysis.IMetricCollector>());
        var manager = new TransientDatabaseResultsManager(
            aggregator,
            Path.Combine(_testOutputDir, "cache.csv"));

        bool removed = manager.RemoveFromCache("NonExistent");

        Assert.That(removed, Is.False);
    }

    private sealed class ThrowingCollector : TaskLayer.ParallelSearch.Analysis.IMetricCollector
    {
        public string CollectorName => "Thrower";

        public IEnumerable<string> GetOutputColumns()
        {
            yield return "Unreachable";
        }

        public Dictionary<string, object> CollectData(TaskLayer.ParallelSearch.Analysis.TransientDatabaseContext context)
        {
            throw new InvalidOperationException("Expected test failure");
        }

        public bool CanCollectData(TaskLayer.ParallelSearch.Analysis.TransientDatabaseContext context)
        {
            return true;
        }
    }
}