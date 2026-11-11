using System;
using System.Collections.Generic;
using NUnit.Framework;
using TaskLayer.ParallelSearch.Analysis;
using TaskLayer.ParallelSearch.Analysis.Collectors;
using Test.ParallelSearchTask.Utility;

namespace Test.ParallelSearchTask.Analysis;

[TestFixture]
public class MetricAggregatorTests
{
    [Test]
    public void RunAnalysis_WithBasicCollector_PopulatesTypedProperties()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var psm = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 100, psmQValue: 0.001, peptideQValue: 0.001);

        var context = ParallelSearchTestContextFactory.CreateContext(
            commonParameters,
            allPsms: [psm],
            transientPsms: [psm],
            allPeptides: [psm],
            transientPeptides: [psm],
            totalProteins: 5,
            transientPeptideCount: 1);

        var aggregator = new MetricAggregator([new BasicMetricCollector()]);

        var result = aggregator.RunAnalysis(context);

        Assert.Multiple(() =>
        {
            Assert.That(result.DatabaseName, Is.EqualTo("TransientDb"));
            Assert.That(result.TargetPsmsAtQValueThreshold, Is.EqualTo(1));
            Assert.That(result.TargetPeptidesAtQValueThreshold, Is.EqualTo(1));
            Assert.That(result.TotalProteins, Is.EqualTo(5));
            Assert.That(result.Errors, Is.Empty);
        });
    }

    [Test]
    public void RunAnalysis_WhenCollectorThrows_RecordsErrorAndContinues()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var psm = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 100, psmQValue: 0.001, peptideQValue: 0.001);

        var context = ParallelSearchTestContextFactory.CreateContext(
            commonParameters,
            allPsms: [psm],
            transientPsms: [psm],
            allPeptides: [psm],
            transientPeptides: [psm]);

        var aggregator = new MetricAggregator([new ThrowingCollector(), new BasicMetricCollector()]);

        var result = aggregator.RunAnalysis(context);

        Assert.Multiple(() =>
        {
            Assert.That(result.Errors.Count, Is.EqualTo(1));
            Assert.That(result.Errors[0], Does.Contain("Thrower"));
            Assert.That(result.TargetPsmsAtQValueThreshold, Is.EqualTo(1));
        });
    }

    private sealed class ThrowingCollector : IMetricCollector
    {
        public string CollectorName => "Thrower";

        public IEnumerable<string> GetOutputColumns()
        {
            yield return "Unreachable";
        }

        public Dictionary<string, object> CollectData(TransientDatabaseContext context)
        {
            throw new InvalidOperationException("Expected test failure");
        }

        public bool CanCollectData(TransientDatabaseContext context)
        {
            return true;
        }
    }
}
