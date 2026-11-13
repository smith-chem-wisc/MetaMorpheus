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
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters(pepQValueThreshold: 1.0); // q-value filtering, as parallel search runs by default

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
            CreateAlignmentBaselineGroup(bestPeptideQValue: 0.05, bestPeptideScore: 200, qValue: 0.01, cumulativeTarget: 100, cumulativeDecoy: 1),
            CreateAlignmentBaselineGroup(bestPeptideQValue: 0.2, bestPeptideScore: 50, qValue: 0.05, cumulativeTarget: 80, cumulativeDecoy: 2),
        };

        var baselineFdr = BaselineFdr.CreateProteinGroupLookup(alignmentBaselineGroups, FilterType.QValue);

        var engine = new TransientProteinScoringAndFdrEngine(
            new List<ProteinGroup> { baselineGroup },
            new List<ProteinGroup> { transientGroup },
            new List<SpectralMatch> { transientPsm },
            baselineFdr,
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
            // its own best peptide, which is what ranked it, is kept
            Assert.That(scoredTransientGroup.BestPeptideQValue, Is.EqualTo(0.12).Within(1e-10));
            Assert.That(scoredTransientGroup.BestPeptideScore, Is.EqualTo(120).Within(1e-10));
            Assert.That(scoredTransientGroup.CumulativeTarget, Is.EqualTo(100));
            Assert.That(scoredTransientGroup.CumulativeDecoy, Is.EqualTo(1));
        });
    }

    [Test]
    public void RunSpecific_WithoutBaselineCache_FallsBackToClassicProteinFdr()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var (transientGroup, transientPsm) = CreateTransientGroupWithPsm(commonParameters, score: 95, psmQValue: 0.04);

        var baselineFdr = BaselineFdr.CreateProteinGroupLookup([], FilterType.QValue);

        var engine = new TransientProteinScoringAndFdrEngine(
            baselineProteinGroups: [],
            transientProteinGroups: [transientGroup],
            neighborhoodPsms: [transientPsm],
            baselineFdr,
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

    /// <summary>
    /// Parsimony hands transient groups over in no particular order. Here the weaker group comes first;
    /// each group must still take the baseline values for its own rank (best peptide q-value, then score).
    /// </summary>
    [Test]
    public void RunSpecific_AlignsEachTransientGroupByItsOwnRank_WhateverOrderParsimonyGaveThem()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters(pepQValueThreshold: 1.0);
        var (weakGroup, weakPsm) = CreateDistinctTransientGroupWithPsm(commonParameters, "PEPTIDEK", "WEAK", score: 20, scanNumber: 601);
        var (strongGroup, strongPsm) = CreateDistinctTransientGroupWithPsm(commonParameters, "ELVISLIVESK", "STRONG", score: 120, scanNumber: 602);

        var baselineFdr = BaselineFdr.CreateProteinGroupLookup(
        [
            CreateAlignmentBaselineGroup(bestPeptideQValue: 0.001, bestPeptideScore: 100, qValue: 0.01, cumulativeTarget: 100, cumulativeDecoy: 1),
            CreateAlignmentBaselineGroup(bestPeptideQValue: 0.001, bestPeptideScore: 30, qValue: 0.05, cumulativeTarget: 80, cumulativeDecoy: 2),
        ], FilterType.QValue);

        var engine = new TransientProteinScoringAndFdrEngine(
            baselineProteinGroups: [],
            transientProteinGroups: [weakGroup, strongGroup],
            neighborhoodPsms: [weakPsm, strongPsm],
            baselineFdr,
            noOneHitWonders: false,
            treatModPeptidesAsDifferentPeptides: true,
            mergeIndistinguishableProteinGroups: false,
            commonParameters,
            fileSpecificParameters: null,
            nestedIds: []);

        var results = (ProteinScoringAndFdrResults)engine.Run();

        Assert.Multiple(() =>
        {
            Assert.That(results.SortedAndScoredProteinGroups, Has.Count.EqualTo(2));
            Assert.That(strongGroup.ProteinGroupScore, Is.GreaterThan(weakGroup.ProteinGroupScore));
            Assert.That(strongGroup.QValue, Is.EqualTo(0.01).Within(1e-10));
            Assert.That(strongGroup.CumulativeTarget, Is.EqualTo(100));
            Assert.That(weakGroup.QValue, Is.EqualTo(0.05).Within(1e-10));
            Assert.That(weakGroup.CumulativeTarget, Is.EqualTo(80));
        });
    }

    /// <summary>
    /// A PEP q-value threshold below the q-value threshold makes protein FDR rank groups by PEP; a lookup built
    /// for q-value ranking must not be used to rank them.
    /// </summary>
    [Test]
    public void RunSpecific_WithALookupRankedByADifferentFilter_Throws()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters(qValueThreshold: 0.01, pepQValueThreshold: 0.005);
        var (group, psm) = CreateDistinctTransientGroupWithPsm(commonParameters, "PEPTIDEK", "ONLY", score: 20, scanNumber: 603);
        var baselineFdr = BaselineFdr.CreateProteinGroupLookup(
            [CreateAlignmentBaselineGroup(bestPeptideQValue: 0.001, bestPeptideScore: 100, qValue: 0.01, cumulativeTarget: 100, cumulativeDecoy: 1)],
            FilterType.QValue);

        var engine = new TransientProteinScoringAndFdrEngine(
            baselineProteinGroups: [],
            transientProteinGroups: [group],
            neighborhoodPsms: [psm],
            baselineFdr,
            noOneHitWonders: false,
            treatModPeptidesAsDifferentPeptides: true,
            mergeIndistinguishableProteinGroups: false,
            commonParameters,
            fileSpecificParameters: null,
            nestedIds: []);

        Assert.Multiple(() =>
        {
            Assert.That(engine.FilterType, Is.EqualTo(FilterType.PepQValue));
            Assert.That(() => engine.Run(), Throws.TypeOf<MetaMorpheusException>());
        });
    }

    /// <summary>
    /// Both constructors pick the filter the same way: PEP q-value only when its threshold is strictly lower,
    /// so equal thresholds rank by q-value, as a regular search does.
    /// </summary>
    [Test]
    [TestCase(0.01, 0.01, FilterType.QValue)]
    [TestCase(0.01, 1.0, FilterType.QValue)]
    [TestCase(0.01, 0.005, FilterType.PepQValue)]
    public void FilterType_IsTheSameWithOrWithoutPreFilteredPsms(double qValueThreshold, double pepQValueThreshold, FilterType expected)
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters(qValueThreshold, pepQValueThreshold);

        var filtering = new ProteinScoringAndFdrEngine([], new List<SpectralMatch>(), false, true, false, commonParameters, null, []);
        var preFiltered = new ProteinScoringAndFdrEngine([], new List<SpectralMatch>(), psmsAlreadyFiltered: true, false, true, false, commonParameters, null, []);

        Assert.Multiple(() =>
        {
            Assert.That(filtering.FilterType, Is.EqualTo(expected));
            Assert.That(preFiltered.FilterType, Is.EqualTo(expected));
        });
    }

    private static (ProteinGroup Group, SpectralMatch Psm) CreateDistinctTransientGroupWithPsm(
        CommonParameters commonParameters, string sequence, string accession, double score, int scanNumber)
    {
        var protein = new Proteomics.Protein(sequence, accession);
        var peptide = protein.Digest(new Proteomics.ProteolyticDigestion.DigestionParams(minPeptideLength: 1), [], [])
            .First(p => p.BaseSequence == sequence);
        var scan = ParallelSearchTestContextFactory.CreateScan(scanNumber, commonParameters);
        var psm = new PeptideSpectralMatch(peptide, 0, score, scanNumber, scan, commonParameters, new List<Omics.Fragmentation.MatchedFragmentIon>());
        psm.ResolveAllAmbiguities();
        psm.SetFdrValues(0, 0, 0.001, 0, 0, 0.001, 0, 0.001);

        var group = new ProteinGroup(
            new HashSet<IBioPolymer> { protein },
            new HashSet<IBioPolymerWithSetMods> { peptide },
            new HashSet<IBioPolymerWithSetMods> { peptide });
        return (group, psm);
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

    private static ProteinGroup CreateAlignmentBaselineGroup(double bestPeptideQValue, double bestPeptideScore, double qValue,
        int cumulativeTarget, int cumulativeDecoy)
    {
        var group = ParallelSearchTestContextFactory.CreateProteinGroup(isDecoy: false, qValue: qValue, peptideCount: 1);
        group.BestPeptideQValue = bestPeptideQValue;
        group.BestPeptideScore = bestPeptideScore;
        group.CumulativeTarget = cumulativeTarget;
        group.CumulativeDecoy = cumulativeDecoy;
        return group;
    }
}
