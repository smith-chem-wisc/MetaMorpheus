using System;
using System.Collections.Generic;
using NUnit.Framework;
using TaskLayer.ParallelSearch.Analysis;
using TaskLayer.ParallelSearch.Analysis.Collectors;
using Test.ParallelSearchTask.Utility;

namespace Test.ParallelSearchTask.Analysis;

[TestFixture]
public class PsmPeptideSearchCollectorTests
{
    [Test]
    public void CanCollectData_WithMissingPeptides_ReturnsFalse()
    {
        var context = new TransientDatabaseContext
        {
            AllPsms = new List<EngineLayer.SpectralMatch>(),
            TransientPsms = new List<EngineLayer.SpectralMatch>(),
            AllPeptides = null!,
            TransientPeptides = new List<EngineLayer.SpectralMatch>(),
            CommonParameters = ParallelSearchTestContextFactory.CreateCommonParameters(),
            TransientDatabase = null!,
            TransientProteins = null!,
            TransientProteinAccessions = new HashSet<string>(),
        };

        var collector = new PsmPeptideSearchCollector();

        Assert.That(collector.CanCollectData(context), Is.False);
    }

    [Test]
    public void CollectData_FiltersByQValueAndReturnsScoreArrays()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters(qValueThreshold: 0.01, pepQValueThreshold: 0.01);

        var targetPass = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 110, psmQValue: 0.005, peptideQValue: 0.005, scanNumber: 1);
        var decoyPass = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: true, score: 90, psmQValue: 0.006, peptideQValue: 0.006, scanNumber: 2);
        var targetFail = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 80, psmQValue: 0.02, peptideQValue: 0.02, scanNumber: 3);

        var allPsms = new List<EngineLayer.SpectralMatch> { targetPass, decoyPass, targetFail };
        var transientPsms = new List<EngineLayer.SpectralMatch> { targetPass, decoyPass };

        var peptideTargetPass = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 105, psmQValue: 0.5, peptideQValue: 0.004, scanNumber: 11);
        var peptideDecoyPass = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: true, score: 70, psmQValue: 0.5, peptideQValue: 0.006, scanNumber: 12);
        var peptideTargetFail = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 60, psmQValue: 0.5, peptideQValue: 0.03, scanNumber: 13);

        var allPeptides = new List<EngineLayer.SpectralMatch> { peptideTargetPass, peptideDecoyPass, peptideTargetFail };
        var transientPeptides = new List<EngineLayer.SpectralMatch> { peptideTargetPass, peptideDecoyPass };

        var context = ParallelSearchTestContextFactory.CreateContext(
            commonParameters,
            allPsms,
            transientPsms,
            allPeptides,
            transientPeptides,
            totalProteins: 2,
            transientPeptideCount: 2);

        var collector = new PsmPeptideSearchCollector();
        var results = collector.CollectData(context);

        Assert.Multiple(() =>
        {
            Assert.That(results[PsmPeptideSearchCollector.PsmTargets], Is.EqualTo(1));
            Assert.That(results[PsmPeptideSearchCollector.PsmDecoys], Is.EqualTo(1));
            Assert.That(results[PsmPeptideSearchCollector.PsmBacterialTargets], Is.EqualTo(1));
            Assert.That(results[PsmPeptideSearchCollector.PsmBacterialDecoys], Is.EqualTo(1));
            Assert.That(results[PsmPeptideSearchCollector.PsmBacterialUnambiguousTargets], Is.EqualTo(1));
            Assert.That(results[PsmPeptideSearchCollector.PsmBacterialUnambiguousDecoys], Is.EqualTo(1));

            Assert.That(results[PsmPeptideSearchCollector.PeptideTargets], Is.EqualTo(1));
            Assert.That(results[PsmPeptideSearchCollector.PeptideDecoys], Is.EqualTo(1));
            Assert.That(results[PsmPeptideSearchCollector.PeptideBacterialTargets], Is.EqualTo(1));
            Assert.That(results[PsmPeptideSearchCollector.PeptideBacterialDecoys], Is.EqualTo(1));
            Assert.That(results[PsmPeptideSearchCollector.PeptideBacterialUnambiguousTargets], Is.EqualTo(1));
            Assert.That(results[PsmPeptideSearchCollector.PeptideBacterialUnambiguousDecoys], Is.EqualTo(1));

            Assert.That((double[])results[PsmPeptideSearchCollector.PsmBacterialUnambiguousTargetScores], Is.EqualTo(new[] { 110d }));
            Assert.That((double[])results[PsmPeptideSearchCollector.PsmBacterialUnambiguousDecoyScores], Is.EqualTo(new[] { 90d }));
            Assert.That((double[])results[PsmPeptideSearchCollector.PeptideBacterialUnambiguousTargetScores], Is.EqualTo(new[] { 105d }));
            Assert.That((double[])results[PsmPeptideSearchCollector.PeptideBacterialUnambiguousDecoyScores], Is.EqualTo(new[] { 70d }));
        });
    }
}
