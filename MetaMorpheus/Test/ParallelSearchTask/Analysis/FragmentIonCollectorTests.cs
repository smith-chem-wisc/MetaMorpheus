using System;
using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using NUnit.Framework;
using TaskLayer.ParallelSearch.Analysis.Collectors;
using Test.ParallelSearchTask.Utility;

namespace Test.ParallelSearchTask.Analysis;

public class FragmentIonCollectorTests
{
    [Test]
    public void CanCollectData_WithMissingPeptides_ReturnsFalse()
    {
        var collector = new FragmentIonCollector();
        var context = ParallelSearchTestContextFactory.CreateContext(
            ParallelSearchTestContextFactory.CreateCommonParameters(),
            allPsms: [],
            transientPsms: [],
            allPeptides: [],
            transientPeptides: null);

        bool result = collector.CanCollectData(context);

        Assert.That(result, Is.False);
    }

    [Test]
    public void CanCollectData_WithMissingPsms_ReturnsFalse()
    {
        var collector = new FragmentIonCollector();
        var context = ParallelSearchTestContextFactory.CreateContext(
            ParallelSearchTestContextFactory.CreateCommonParameters(),
            allPsms: [],
            transientPsms: null,
            allPeptides: [],
            transientPeptides: []);

        bool result = collector.CanCollectData(context);

        Assert.That(result, Is.False);
    }

    [Test]
    public void CanCollectData_WithValidData_ReturnsTrue()
    {
        var collector = new FragmentIonCollector();
        var psm = ParallelSearchTestContextFactory.CreateSpectralMatch(
            ParallelSearchTestContextFactory.CreateCommonParameters(),
            isDecoy: false,
            score: 5.0,
            psmQValue: 0.001,
            peptideQValue: 0.001);

        var context = ParallelSearchTestContextFactory.CreateContext(
            ParallelSearchTestContextFactory.CreateCommonParameters(),
            allPsms: [psm],
            transientPsms: [psm],
            allPeptides: [psm],
            transientPeptides: [psm]);

        bool result = collector.CanCollectData(context);

        Assert.That(result, Is.True);
    }

    [Test]
    public void CanCollectData_WithEmptyLists_ReturnsTrue()
    {
        var collector = new FragmentIonCollector();
        var context = ParallelSearchTestContextFactory.CreateContext(
            ParallelSearchTestContextFactory.CreateCommonParameters(),
            allPsms: [],
            transientPsms: [],
            allPeptides: [],
            transientPeptides: []);

        bool result = collector.CanCollectData(context);

        Assert.That(result, Is.True);
    }

    [Test]
    public void CollectData_WithNoConfidentEntries_ReturnsNaNMedians()
    {
        var collector = new FragmentIonCollector();
        var psmHighQ = ParallelSearchTestContextFactory.CreateSpectralMatch(
            ParallelSearchTestContextFactory.CreateCommonParameters(),
            isDecoy: false,
            score: 5.0,
            psmQValue: 0.5,
            peptideQValue: 0.5);

        var context = ParallelSearchTestContextFactory.CreateContext(
            ParallelSearchTestContextFactory.CreateCommonParameters(qValueThreshold: 0.01),
            allPsms: [psmHighQ],
            transientPsms: [psmHighQ],
            allPeptides: [psmHighQ],
            transientPeptides: [psmHighQ]);

        var result = collector.CollectData(context);

        Assert.Multiple(() =>
        {
            Assert.That(result[FragmentIonCollector.PSM_LongestIonSeriesBidirectionalTargets], Is.EqualTo(double.NaN));
            Assert.That(result[FragmentIonCollector.PSM_ComplementaryIonCountTargets], Is.EqualTo(double.NaN));
            Assert.That(result[FragmentIonCollector.PSM_SequenceCoverageFractionTargets], Is.EqualTo(double.NaN));
        });
    }

    [Test]
    public void CollectData_FiltersByQValueAndReturnsMetrics()
    {
        var commonParams = ParallelSearchTestContextFactory.CreateCommonParameters(qValueThreshold: 0.05);

        var targetPsm = ParallelSearchTestContextFactory.CreateSpectralMatch(
            commonParams,
            isDecoy: false,
            score: 10.0,
            psmQValue: 0.001,
            peptideQValue: 0.001);

        var decoyPsm = ParallelSearchTestContextFactory.CreateSpectralMatch(
            commonParams,
            isDecoy: true,
            score: 8.0,
            psmQValue: 0.001,
            peptideQValue: 0.001);

        var context = ParallelSearchTestContextFactory.CreateContext(
            commonParams,
            allPsms: [targetPsm, decoyPsm],
            transientPsms: [targetPsm, decoyPsm],
            allPeptides: [targetPsm, decoyPsm],
            transientPeptides: [targetPsm, decoyPsm]);

        var collector = new FragmentIonCollector();
        var result = collector.CollectData(context);

        Assert.Multiple(() =>
        {
            Assert.That(result[FragmentIonCollector.PSM_LongestIonSeriesBidirectionalTargets], Is.Not.Null);
            Assert.That(result[FragmentIonCollector.PSM_LongestIonSeriesBidirectionalDecoys], Is.Not.Null);
            Assert.That(result[FragmentIonCollector.PSM_LongestIonSeriesBidirectional_AllTargets], Is.TypeOf<double[]>());
        });
    }

    [Test]
    public void CollectData_ReturnsAllColumnNames()
    {
        var collector = new FragmentIonCollector();
        var columns = collector.GetOutputColumns().ToList();

        Assert.Multiple(() =>
        {
            Assert.That(columns, Contains.Item(FragmentIonCollector.PSM_LongestIonSeriesBidirectionalTargets));
            Assert.That(columns, Contains.Item(FragmentIonCollector.PSM_ComplementaryIonCountTargets));
            Assert.That(columns, Contains.Item(FragmentIonCollector.PSM_SequenceCoverageFractionTargets));
            Assert.That(columns, Contains.Item(FragmentIonCollector.PSM_LongestIonSeriesBidirectionalDecoys));
            Assert.That(columns, Contains.Item(FragmentIonCollector.PSM_ComplementaryIonCountDecoys));
            Assert.That(columns, Contains.Item(FragmentIonCollector.PSM_SequenceCoverageFractionDecoys));
            Assert.That(columns, Contains.Item(FragmentIonCollector.Peptide_LongestIonSeriesBidirectionalTargets));
            Assert.That(columns, Contains.Item(FragmentIonCollector.Peptide_LongestIonSeriesBidirectionalDecoys));
        });

        Assert.That(columns.Count, Is.EqualTo(24));
    }
}