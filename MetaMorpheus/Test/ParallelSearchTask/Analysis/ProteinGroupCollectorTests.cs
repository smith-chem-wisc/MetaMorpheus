using System.Collections.Generic;
using NUnit.Framework;
using TaskLayer.ParallelSearch.Analysis;
using TaskLayer.ParallelSearch.Analysis.Collectors;
using Test.ParallelSearchTask.Utility;

namespace Test.ParallelSearchTask.Analysis;

[TestFixture]
public class ProteinGroupCollectorTests
{
    [Test]
    public void CanCollectData_WhenProteinGroupsMissing_ReturnsFalse()
    {
        var context = ParallelSearchTestContextFactory.CreateContext(
            ParallelSearchTestContextFactory.CreateCommonParameters(),
            [],
            [],
            [],
            []);

        var collector = new ProteinGroupCollector();

        Assert.That(collector.CanCollectData(context), Is.False);
    }

    [Test]
    public void CollectData_FiltersByQValueAndMinimumPeptides()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters(qValueThreshold: 0.01, pepQValueThreshold: 0.05);

        var passingTarget = ParallelSearchTestContextFactory.CreateProteinGroup(isDecoy: false, qValue: 0.005, peptideCount: 2);
        var passingDecoy = ParallelSearchTestContextFactory.CreateProteinGroup(isDecoy: true, qValue: 0.008, peptideCount: 2);
        var failingByQValue = ParallelSearchTestContextFactory.CreateProteinGroup(isDecoy: false, qValue: 0.03, peptideCount: 2);
        var failingByPeptideCount = ParallelSearchTestContextFactory.CreateProteinGroup(isDecoy: false, qValue: 0.001, peptideCount: 1);

        var context = ParallelSearchTestContextFactory.CreateContext(
            commonParameters,
            [],
            [],
            [],
            [],
            proteinGroups: [passingTarget, passingDecoy, failingByQValue, failingByPeptideCount],
            transientProteinGroups: [passingTarget, passingDecoy, failingByQValue, failingByPeptideCount]);

        var collector = new ProteinGroupCollector();
        var results = collector.CollectData(context);

        Assert.Multiple(() =>
        {
            Assert.That(results[ProteinGroupCollector.ProteinGroupTargets], Is.EqualTo(3));
            Assert.That(results[ProteinGroupCollector.ProteinGroupDecoys], Is.EqualTo(1));
            Assert.That(results[ProteinGroupCollector.TargetProteinGroupsAtQValueThreshold], Is.EqualTo(1));
            Assert.That(results[ProteinGroupCollector.TargetProteinGroupsFromTransientDb], Is.EqualTo(3));
            Assert.That(results[ProteinGroupCollector.TargetProteinGroupsFromTransientDbAtQValueThreshold], Is.EqualTo(1));
            Assert.That(results[ProteinGroupCollector.ProteinGroupBacterialTargets], Is.EqualTo(3));
            Assert.That(results[ProteinGroupCollector.ProteinGroupBacterialDecoys], Is.EqualTo(1));
            Assert.That(results[ProteinGroupCollector.ProteinGroupBacterialUnambiguousTargets], Is.EqualTo(1));
            Assert.That(results[ProteinGroupCollector.ProteinGroupBacterialUnambiguousDecoys], Is.EqualTo(1));
        });
    }
}
