using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using EngineLayer.ParallelSearch;
using EngineLayer.ParallelSearch.FdrAlignment;
using EngineLayer.SpectrumMatch;
using NUnit.Framework;
using Omics;
using Test.ParallelSearchTask.Utility;

namespace Test.ParallelSearchTask.Engine;

[TestFixture]
public class TransientProteinScoringAndFdrEngineTests
{
    [Test]
    public void RunSpecific_UsesBaselineAlignmentForTransientGroups()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();

        var baselineGroup = ParallelSearchTestContextFactory.CreateProteinGroup(isDecoy: false, qValue: 0.33, peptideCount: 1);
        baselineGroup.ProteinGroupScore = 250;
        baselineGroup.BestPeptideScore = 250;
        baselineGroup.BestPeptideQValue = 0.33;
        baselineGroup.BestPeptidePEP = 0.33;
        baselineGroup.CumulativeTarget = 11;
        baselineGroup.CumulativeDecoy = 2;

        var (transientGroup, transientPsm) = CreateTransientGroupWithPsm(commonParameters, score: 120, psmQValue: 0.12);

        var alignmentBaselineGroups = new List<ProteinGroup>
        {
            CreateAlignmentBaselineGroup(score: 200, qValue: 0.01, pepQValue: 0.02, pep: 0.03, cumulativeTarget: 100, cumulativeDecoy: 1),
            CreateAlignmentBaselineGroup(score: 50, qValue: 0.05, pepQValue: 0.06, pep: 0.07, cumulativeTarget: 80, cumulativeDecoy: 2),
        };

        var alignmentService = new ProteinGroupFdrAlignmentService();
        alignmentService.BuildBaselineCache(alignmentBaselineGroups);

        var engine = new TransientProteinScoringAndFdrEngine(
            new List<ProteinGroup> { baselineGroup },
            new List<ProteinGroup> { transientGroup },
            new List<SpectralMatch> { transientPsm },
            alignmentService,
            noOneHitWonders: false,
            treatModPeptidesAsDifferentPeptides: true,
            mergeIndistinguishableProteinGroups: true,
            commonParameters,
            fileSpecificParameters: null,
            nestedIds: []);

        var results = (ProteinScoringAndFdrResults)engine.Run();
        var scoredTransientGroup = results.SortedAndScoredProteinGroups.Single(pg => pg != baselineGroup);

        Assert.Multiple(() =>
        {
            Assert.That(results.SortedAndScoredProteinGroups.Count, Is.EqualTo(2));
            Assert.That(baselineGroup.QValue, Is.EqualTo(0.33).Within(1e-10));
            Assert.That(scoredTransientGroup.ProteinGroupScore, Is.EqualTo(120).Within(1e-10));
            Assert.That(scoredTransientGroup.AllPsmsBelowOnePercentFDR.Count, Is.EqualTo(1));
            Assert.That(scoredTransientGroup.QValue, Is.EqualTo(0.01).Within(1e-10));
            Assert.That(scoredTransientGroup.BestPeptideQValue, Is.EqualTo(0.02).Within(1e-10));
            Assert.That(scoredTransientGroup.BestPeptidePEP, Is.EqualTo(0.03).Within(1e-10));
            Assert.That(scoredTransientGroup.CumulativeTarget, Is.EqualTo(100));
            Assert.That(scoredTransientGroup.CumulativeDecoy, Is.EqualTo(1));
        });
    }

    [Test]
    public void RunSpecific_WithoutBaselineCache_FallsBackToClassicProteinFdr()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var (transientGroup, transientPsm) = CreateTransientGroupWithPsm(commonParameters, score: 95, psmQValue: 0.04);

        var alignmentService = new ProteinGroupFdrAlignmentService();
        alignmentService.BuildBaselineCache([]);

        var engine = new TransientProteinScoringAndFdrEngine(
            baselineProteinGroups: [],
            transientProteinGroups: [transientGroup],
            neighborhoodPsms: [transientPsm],
            alignmentService,
            noOneHitWonders: false,
            treatModPeptidesAsDifferentPeptides: true,
            mergeIndistinguishableProteinGroups: true,
            commonParameters,
            fileSpecificParameters: null,
            nestedIds: []);

        var results = (ProteinScoringAndFdrResults)engine.Run();
        var scoredTransientGroup = results.SortedAndScoredProteinGroups.Single();

        Assert.Multiple(() =>
        {
            Assert.That(results.SortedAndScoredProteinGroups.Count, Is.EqualTo(1));
            Assert.That(scoredTransientGroup.ProteinGroupScore, Is.EqualTo(95).Within(1e-10));
            Assert.That(scoredTransientGroup.CumulativeTarget, Is.EqualTo(1));
            Assert.That(scoredTransientGroup.CumulativeDecoy, Is.EqualTo(0));
            Assert.That(scoredTransientGroup.QValue, Is.EqualTo(0).Within(1e-10));
            Assert.That(scoredTransientGroup.BestPeptideQValue, Is.EqualTo(0.04).Within(1e-10));
            Assert.That(scoredTransientGroup.BestPeptidePEP, Is.EqualTo(0).Within(1e-10));
        });
    }

    private static (ProteinGroup Group, SpectralMatch Psm) CreateTransientGroupWithPsm(
        CommonParameters commonParameters,
        double score,
        double psmQValue)
    {
        var psm = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, false, score, psmQValue, psmQValue, 501);
        var peptide = psm.BestMatchingBioPolymersWithSetMods.Single().SpecificBioPolymer;
        var transientGroup = new ProteinGroup(
            new HashSet<IBioPolymer> { peptide.Parent },
            new HashSet<IBioPolymerWithSetMods> { peptide },
            new HashSet<IBioPolymerWithSetMods> { peptide });

        return (transientGroup, psm);
    }

    private static ProteinGroup CreateAlignmentBaselineGroup(double score, double qValue, double pepQValue, double pep,
        int cumulativeTarget, int cumulativeDecoy)
    {
        var group = ParallelSearchTestContextFactory.CreateProteinGroup(isDecoy: false, qValue: qValue, peptideCount: 1);
        group.ProteinGroupScore = score;
        group.BestPeptideQValue = pepQValue;
        group.BestPeptidePEP = pep;
        group.CumulativeTarget = cumulativeTarget;
        group.CumulativeDecoy = cumulativeDecoy;
        return group;
    }
}
