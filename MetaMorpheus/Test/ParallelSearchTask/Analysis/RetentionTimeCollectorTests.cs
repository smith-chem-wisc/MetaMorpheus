using System;
using System.Collections.Generic;
using NUnit.Framework;
using TaskLayer.ParallelSearch.Analysis;
using TaskLayer.ParallelSearch.Analysis.Collectors;
using Test.ParallelSearchTask.Utility;

namespace Test.ParallelSearchTask.Analysis;

[TestFixture]
public class RetentionTimeCollectorTests
{
    [Test]
    public void CanCollectData_WithMissingLists_ReturnsFalse()
    {
        var context = new TransientDatabaseContext
        {
            AllPsms = null!,
            TransientPsms = new List<EngineLayer.SpectralMatch>(),
            AllPeptides = new List<EngineLayer.SpectralMatch>(),
            TransientPeptides = new List<EngineLayer.SpectralMatch>(),
            CommonParameters = ParallelSearchTestContextFactory.CreateCommonParameters(),
        };

        var collector = new RetentionTimeCollector();

        Assert.That(collector.CanCollectData(context), Is.False);
    }

    [Test]
    public void CollectData_WithNoConfidentEntries_ReturnsNaNAndEmptyArrays()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters(qValueThreshold: 0.01, pepQValueThreshold: 0.01);
        var highQValuePsm = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 100, psmQValue: 0.5, peptideQValue: 0.5);

        var context = ParallelSearchTestContextFactory.CreateContext(
            commonParameters,
            allPsms: [highQValuePsm],
            transientPsms: [highQValuePsm],
            allPeptides: [highQValuePsm],
            transientPeptides: [highQValuePsm]);

        var collector = new RetentionTimeCollector();
        var result = collector.CollectData(context);

        Assert.Multiple(() =>
        {
            Assert.That(double.IsNaN((double)result[RetentionTimeCollector.PsmMeanAbsoluteRtError]), Is.True);
            Assert.That(double.IsNaN((double)result[RetentionTimeCollector.PsmRtCorrelationCoefficient]), Is.True);
            Assert.That(double.IsNaN((double)result[RetentionTimeCollector.PeptideMeanAbsoluteRtError]), Is.True);
            Assert.That(double.IsNaN((double)result[RetentionTimeCollector.PeptideRtCorrelationCoefficient]), Is.True);
            Assert.That((double[])result[RetentionTimeCollector.PsmAllRtErrors], Is.Empty);
            Assert.That((double[])result[RetentionTimeCollector.PeptideAllRtErrors], Is.Empty);
        });
    }

    [Test]
    public void CollectData_FiltersPsmsByPsmQValueAndPeptidesByPeptideQValue()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters(qValueThreshold: 0.01, pepQValueThreshold: 0.01);

        // passes at the PSM level only
        var psmLevelOnly = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 100, psmQValue: 0.001, peptideQValue: 0.5, scanNumber: 1);

        var context = ParallelSearchTestContextFactory.CreateContext(
            commonParameters,
            allPsms: [psmLevelOnly],
            transientPsms: [psmLevelOnly],
            allPeptides: [psmLevelOnly],
            transientPeptides: [psmLevelOnly]);

        var result = new RetentionTimeCollector().CollectData(context);

        Assert.Multiple(() =>
        {
            Assert.That((double[])result[RetentionTimeCollector.PsmAllRtErrors], Has.Length.EqualTo(1));
            Assert.That((double[])result[RetentionTimeCollector.PeptideAllRtErrors], Is.Empty);
            Assert.That(double.IsNaN((double)result[RetentionTimeCollector.PeptideMeanAbsoluteRtError]), Is.True);
        });

        // passes at the peptide level only
        var peptideLevelOnly = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 100, psmQValue: 0.5, peptideQValue: 0.001, scanNumber: 2);

        context = ParallelSearchTestContextFactory.CreateContext(
            commonParameters,
            allPsms: [peptideLevelOnly],
            transientPsms: [peptideLevelOnly],
            allPeptides: [peptideLevelOnly],
            transientPeptides: [peptideLevelOnly]);

        result = new RetentionTimeCollector().CollectData(context);

        Assert.Multiple(() =>
        {
            Assert.That((double[])result[RetentionTimeCollector.PsmAllRtErrors], Is.Empty);
            Assert.That((double[])result[RetentionTimeCollector.PeptideAllRtErrors], Has.Length.EqualTo(1));
        });
    }
}
