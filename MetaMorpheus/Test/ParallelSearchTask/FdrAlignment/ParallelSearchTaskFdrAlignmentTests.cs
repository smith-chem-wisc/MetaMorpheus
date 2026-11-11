using System;
using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.ParallelSearch.FdrAlignment;
using EngineLayer.SpectrumMatch;
using NUnit.Framework;
using Test.ParallelSearchTask.Utility;

namespace Test.ParallelSearchTask.FdrAlignment;

[TestFixture]
public class ParallelSearchTaskFdrAlignmentTests
{
    [Test]
    public void Find_TakesTheLowestBaselineEntryScoringAtLeastAsHigh_AndClampsAtBothEnds()
    {
        var lookup = new BaselineScoreLookup<string>([(100, "a"), (50, "b"), (10, "c")]);

        Assert.Multiple(() =>
        {
            Assert.That(lookup.Count, Is.EqualTo(3));
            Assert.That(lookup.Find(110), Is.EqualTo("a"), "above every baseline score");
            Assert.That(lookup.Find(100), Is.EqualTo("a"));
            Assert.That(lookup.Find(80), Is.EqualTo("a"));
            Assert.That(lookup.Find(50), Is.EqualTo("b"));
            Assert.That(lookup.Find(40), Is.EqualTo("b"));
            Assert.That(lookup.Find(10), Is.EqualTo("c"));
            Assert.That(lookup.Find(5), Is.EqualTo("c"), "below every baseline score");
        });
    }

    [Test]
    public void Find_SortsTheBaselineByScore_AndKeepsTiesInTheOrderGiven()
    {
        var lookup = new BaselineScoreLookup<string>([(10, "c"), (50, "first 50"), (100, "a"), (50, "second 50")]);

        Assert.Multiple(() =>
        {
            Assert.That(lookup.Find(100), Is.EqualTo("a"));
            // every entry scoring >= 50 counts, so the last of the tied entries is taken
            Assert.That(lookup.Find(50), Is.EqualTo("second 50"));
            Assert.That(lookup.Find(60), Is.EqualTo("a"));
            Assert.That(lookup.Find(20), Is.EqualTo("second 50"));
        });
    }

    [Test]
    public void EmptyLookup_LeavesItemsAsTheyAre_AndFindThrows()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var psm = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, false, 40, 0.9, 0.9, 1);
        var originalFdr = psm.PsmFdrInfo;
        var group = CreateTransientGroup(bestPeptideQValue: 0.001, bestPeptideScore: 20);

        var originalPeptideFdr = psm.PeptideFdrInfo;

        var psmLookup = BaselineFdr.CreatePsmLookup([]);
        BaselineFdr.AlignPsms(psmLookup, [psm]);
        BaselineFdr.AlignPeptides(BaselineFdr.CreatePeptideLookup([]), [psm]);
        BaselineFdr.AlignProteinGroups(BaselineFdr.CreateProteinGroupLookup([], FilterType.QValue), [group]);

        Assert.Multiple(() =>
        {
            Assert.That(psmLookup.IsEmpty, Is.True);
            Assert.That(psm.PsmFdrInfo, Is.SameAs(originalFdr));
            Assert.That(psm.PeptideFdrInfo, Is.SameAs(originalPeptideFdr));
            Assert.That(group.QValue, Is.EqualTo(0.9).Within(1e-10));
            Assert.That(() => psmLookup.Find(40), Throws.InvalidOperationException);
        });
    }

    [Test]
    public void AlignPsms_CopiesTheBaselineFdrOfTheMatchingScore()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var baselinePsms = CreateBaselinePsms(commonParameters);
        var lookup = BaselineFdr.CreatePsmLookup(baselinePsms);

        var transientPsms = CreateTransientPsms(commonParameters, 110, 80, 40, 5);
        BaselineFdr.AlignPsms(lookup, transientPsms);

        Assert.Multiple(() =>
        {
            Assert.That(transientPsms[0].PsmFdrInfo.QValue, Is.EqualTo(0.01).Within(1e-10));
            Assert.That(transientPsms[1].PsmFdrInfo.QValue, Is.EqualTo(0.01).Within(1e-10));
            Assert.That(transientPsms[2].PsmFdrInfo.QValue, Is.EqualTo(0.05).Within(1e-10));
            Assert.That(transientPsms[3].PsmFdrInfo.QValue, Is.EqualTo(0.10).Within(1e-10));

            // each transient PSM gets its own copy, so changing one leaves the baseline alone
            Assert.That(transientPsms[0].PsmFdrInfo, Is.Not.SameAs(transientPsms[1].PsmFdrInfo));
            Assert.That(transientPsms[0].PsmFdrInfo, Is.Not.SameAs(baselinePsms[0].PsmFdrInfo));
        });

        transientPsms[0].PsmFdrInfo.QValue = 0.5;
        Assert.That(lookup.Find(110).QValue, Is.EqualTo(0.01).Within(1e-10));
    }

    [Test]
    public void AlignPeptides_UsesOnlyBaselinePsmsWithPeptideFdr()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var baselinePsms = new List<SpectralMatch>
        {
            CreateSpectralMatchWithFdr(commonParameters, 100, psmQValue: 0.01, peptideQValue: null, scanIndex: 111),
            CreateSpectralMatchWithFdr(commonParameters, 50, psmQValue: 0.05, peptideQValue: 0.07, scanIndex: 112),
            CreateSpectralMatchWithFdr(commonParameters, 10, psmQValue: 0.10, peptideQValue: 0.2, scanIndex: 113),
            CreateSpectralMatchWithFdr(commonParameters, 5, psmQValue: 0.10, peptideQValue: 2, scanIndex: 114),
        };
        var lookup = BaselineFdr.CreatePeptideLookup(baselinePsms);

        var transientPeptides = CreateTransientPsms(commonParameters, 100, 40, 1);
        BaselineFdr.AlignPeptides(lookup, transientPeptides);

        Assert.Multiple(() =>
        {
            Assert.That(lookup.Count, Is.EqualTo(2), "no peptide FDR, or a q-value of 2, is left out");
            Assert.That(transientPeptides[0].PeptideFdrInfo.QValue, Is.EqualTo(0.07).Within(1e-10));
            Assert.That(transientPeptides[1].PeptideFdrInfo.QValue, Is.EqualTo(0.07).Within(1e-10));
            Assert.That(transientPeptides[2].PeptideFdrInfo.QValue, Is.EqualTo(0.2).Within(1e-10));
        });
    }

    [Test]
    public void AlignProteinGroups_RanksByBestPeptide_AndCopiesOnlyTheProteinFdr()
    {
        var lookup = BaselineFdr.CreateProteinGroupLookup(CreateBaselineGroups(), FilterType.QValue);

        var transientGroups = new List<ProteinGroup>
        {
            CreateTransientGroup(bestPeptideQValue: 0.0005, bestPeptideScore: 25), // above every baseline group
            CreateTransientGroup(bestPeptideQValue: 0.001, bestPeptideScore: 15),  // ties A on q-value, lower score
            CreateTransientGroup(bestPeptideQValue: 0.004, bestPeptideScore: 12),  // ties B exactly
            CreateTransientGroup(bestPeptideQValue: 0.01, bestPeptideScore: 9),
            CreateTransientGroup(bestPeptideQValue: 0.5, bestPeptideScore: 1),     // below every baseline group
        };
        BaselineFdr.AlignProteinGroups(lookup, transientGroups);

        Assert.Multiple(() =>
        {
            Assert.That(transientGroups.Select(g => g.QValue), Is.EqualTo(new[] { 0.01, 0.01, 0.05, 0.05, 0.10 }).Within(1e-10));
            Assert.That(transientGroups.Select(g => g.CumulativeTarget), Is.EqualTo(new[] { 100, 100, 80, 80, 60 }));
            Assert.That(transientGroups.Select(g => g.CumulativeDecoy), Is.EqualTo(new[] { 1, 1, 2, 2, 3 }));

            // the groups keep their own best-peptide values; those are what placed each group
            Assert.That(transientGroups.Select(g => g.BestPeptideQValue), Is.EqualTo(new[] { 0.0005, 0.001, 0.004, 0.01, 0.5 }).Within(1e-10));
            Assert.That(transientGroups.Select(g => g.BestPeptideScore), Is.EqualTo(new[] { 25.0, 15, 12, 9, 1 }).Within(1e-10));
            Assert.That(transientGroups.Select(g => g.BestPeptidePEP), Has.All.EqualTo(0.9).Within(1e-10));
        });
    }

    /// <summary>
    /// Protein FDR ranks groups by their best peptide, not by the summed protein-group score. Baseline group B
    /// has the highest summed score (many weak peptides) but ranks second. A transient group just below B's
    /// summed score but with a best peptide as good as A's must take A's q-value; keyed on summed score it took B's.
    /// </summary>
    [Test]
    public void AlignProteinGroups_IgnoresTheSummedProteinGroupScore()
    {
        var lookup = BaselineFdr.CreateProteinGroupLookup(CreateBaselineGroups(), FilterType.QValue);
        var transient = CreateTransientGroup(bestPeptideQValue: 0.001, bestPeptideScore: 20.5, summedScore: 199);

        BaselineFdr.AlignProteinGroups(lookup, [transient]);

        Assert.That(transient.QValue, Is.EqualTo(0.01).Within(1e-10));
    }

    [Test]
    public void AlignProteinGroups_WithPepFiltering_RanksByBestPeptidePep()
    {
        var lookup = BaselineFdr.CreateProteinGroupLookup(CreateBaselineGroups(), FilterType.PepQValue);

        // a best-peptide q-value that would place it last, but a PEP that places it first
        var transient = CreateTransientGroup(bestPeptideQValue: 0.5, bestPeptideScore: 20, bestPeptidePep: 0.005);
        BaselineFdr.AlignProteinGroups(lookup, [transient]);

        Assert.Multiple(() =>
        {
            Assert.That(lookup.FilterType, Is.EqualTo(FilterType.PepQValue));
            Assert.That(transient.QValue, Is.EqualTo(0.01).Within(1e-10));
        });
    }

    /// <summary>
    /// Transient protein groups arrive in parsimony order, not rank order, and PSMs need not arrive sorted
    /// either. Each item must get the baseline values for its own rank.
    /// </summary>
    [Test]
    public void Alignment_DoesNotDependOnTheOrderItemsArriveIn()
    {
        var groupLookup = BaselineFdr.CreateProteinGroupLookup(CreateBaselineGroups(), FilterType.QValue);
        var groups = new List<ProteinGroup>
        {
            CreateTransientGroup(bestPeptideQValue: 0.5, bestPeptideScore: 1),
            CreateTransientGroup(bestPeptideQValue: 0.0005, bestPeptideScore: 25),
            CreateTransientGroup(bestPeptideQValue: 0.01, bestPeptideScore: 9),
            CreateTransientGroup(bestPeptideQValue: 0.001, bestPeptideScore: 15),
        };
        BaselineFdr.AlignProteinGroups(groupLookup, groups);

        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var psmLookup = BaselineFdr.CreatePsmLookup(CreateBaselinePsms(commonParameters));
        var psms = CreateTransientPsms(commonParameters, 5, 110, 40, 80);
        BaselineFdr.AlignPsms(psmLookup, psms);

        Assert.Multiple(() =>
        {
            Assert.That(groups.Select(g => g.QValue), Is.EqualTo(new[] { 0.10, 0.01, 0.05, 0.01 }).Within(1e-10));
            Assert.That(psms.Select(p => p.PsmFdrInfo.QValue), Is.EqualTo(new[] { 0.10, 0.01, 0.05, 0.01 }).Within(1e-10));
        });
    }

    private static List<SpectralMatch> CreateBaselinePsms(CommonParameters commonParameters)
    {
        return
        [
            CreateSpectralMatchWithFdr(commonParameters, 100, psmQValue: 0.01, peptideQValue: 0.01, scanIndex: 101),
            CreateSpectralMatchWithFdr(commonParameters, 50, psmQValue: 0.05, peptideQValue: 0.05, scanIndex: 102),
            CreateSpectralMatchWithFdr(commonParameters, 10, psmQValue: 0.10, peptideQValue: 0.10, scanIndex: 103),
        ];
    }

    private static List<SpectralMatch> CreateTransientPsms(CommonParameters commonParameters, params double[] scores)
    {
        return scores
            .Select((score, i) => ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, false, score, 0.9, 0.9, i + 1))
            .ToList();
    }

    /// <summary>
    /// Ranked by best peptide: A, then B, then C. B has by far the highest summed score.
    /// </summary>
    private static List<ProteinGroup> CreateBaselineGroups()
    {
        return
        [
            CreateProteinGroup(bestPeptideQValue: 0.001, bestPeptideScore: 20, bestPeptidePep: 0.01, qValue: 0.01, cumulativeTarget: 100, cumulativeDecoy: 1, summedScore: 30),
            CreateProteinGroup(bestPeptideQValue: 0.004, bestPeptideScore: 12, bestPeptidePep: 0.2, qValue: 0.05, cumulativeTarget: 80, cumulativeDecoy: 2, summedScore: 200),
            CreateProteinGroup(bestPeptideQValue: 0.02, bestPeptideScore: 8, bestPeptidePep: 0.4, qValue: 0.10, cumulativeTarget: 60, cumulativeDecoy: 3, summedScore: 15),
        ];
    }

    private static ProteinGroup CreateTransientGroup(double bestPeptideQValue, double bestPeptideScore,
        double bestPeptidePep = 0.9, double summedScore = 0)
    {
        return CreateProteinGroup(bestPeptideQValue, bestPeptideScore, bestPeptidePep, qValue: 0.9,
            cumulativeTarget: 1, cumulativeDecoy: 1, summedScore);
    }

    private static FdrInfo CreateFdrInfo(double qValue)
    {
        return new FdrInfo
        {
            QValue = qValue,
            QValueNotch = qValue,
            PEP = qValue,
            PEP_QValue = qValue,
            CumulativeDecoy = 1,
            CumulativeTarget = 10,
            CumulativeDecoyNotch = 1,
            CumulativeTargetNotch = 10,
        };
    }

    private static SpectralMatch CreateSpectralMatchWithFdr(CommonParameters commonParameters, double score,
        double psmQValue, double? peptideQValue, int scanIndex)
    {
        var psm = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, false, score, 0.9, 0.9, scanIndex);
        psm.PsmFdrInfo = CreateFdrInfo(psmQValue);
        psm.PeptideFdrInfo = peptideQValue.HasValue ? CreateFdrInfo(peptideQValue.Value) : null;
        return psm;
    }

    private static ProteinGroup CreateProteinGroup(double bestPeptideQValue, double bestPeptideScore, double bestPeptidePep,
        double qValue, int cumulativeTarget, int cumulativeDecoy, double summedScore)
    {
        var proteinGroup = ParallelSearchTestContextFactory.CreateProteinGroup(isDecoy: false, qValue: qValue, peptideCount: 1);
        proteinGroup.ProteinGroupScore = summedScore;
        proteinGroup.BestPeptideQValue = bestPeptideQValue;
        proteinGroup.BestPeptideScore = bestPeptideScore;
        proteinGroup.BestPeptidePEP = bestPeptidePep;
        proteinGroup.CumulativeTarget = cumulativeTarget;
        proteinGroup.CumulativeDecoy = cumulativeDecoy;
        return proteinGroup;
    }
}
