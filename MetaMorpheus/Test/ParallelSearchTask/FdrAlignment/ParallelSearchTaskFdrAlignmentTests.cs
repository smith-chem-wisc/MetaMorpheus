using System;
using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.ParallelSearch.FdrAlignment;
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
        var group = CreateProteinGroupWithFdr(score: 40, qValue: 0.9, bestPeptideQValue: 0.9, bestPeptidePep: 0.9, cumulativeTarget: 1, cumulativeDecoy: 1);

        var originalPeptideFdr = psm.PeptideFdrInfo;

        var psmLookup = BaselineFdr.CreatePsmLookup([]);
        BaselineFdr.AlignPsms(psmLookup, [psm]);
        BaselineFdr.AlignPeptides(BaselineFdr.CreatePeptideLookup([]), [psm]);
        BaselineFdr.AlignProteinGroups(BaselineFdr.CreateProteinGroupLookup([]), [group]);

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
    public void AlignProteinGroups_CopiesEveryBaselineFdrField()
    {
        var lookup = BaselineFdr.CreateProteinGroupLookup(CreateBaselineGroups());

        var transientGroups = new[] { 110.0, 80, 40, 5 }.Select(CreateTransientGroup).ToList();
        BaselineFdr.AlignProteinGroups(lookup, transientGroups);

        Assert.Multiple(() =>
        {
            Assert.That(transientGroups[0].QValue, Is.EqualTo(0.01).Within(1e-10));
            Assert.That(transientGroups[1].QValue, Is.EqualTo(0.01).Within(1e-10));
            Assert.That(transientGroups[2].QValue, Is.EqualTo(0.05).Within(1e-10));
            Assert.That(transientGroups[3].QValue, Is.EqualTo(0.10).Within(1e-10));
            Assert.That(transientGroups[2].BestPeptideQValue, Is.EqualTo(0.06).Within(1e-10));
            Assert.That(transientGroups[3].BestPeptidePEP, Is.EqualTo(0.12).Within(1e-10));
            Assert.That(transientGroups[1].CumulativeTarget, Is.EqualTo(100));
            Assert.That(transientGroups[3].CumulativeDecoy, Is.EqualTo(3));
        });
    }

    /// <summary>
    /// Transient protein groups arrive in parsimony order, not score order. Each one must get the
    /// baseline values for its own score; the old forward-only walk gave a high-scoring group that
    /// came after a low-scoring one the low-scoring group's values.
    /// </summary>
    [Test]
    public void Alignment_DoesNotDependOnTheOrderItemsArriveIn()
    {
        double[] scores = [5, 110, 40, 80];
        double[] expectedQValues = [0.10, 0.01, 0.05, 0.01];

        var groupLookup = BaselineFdr.CreateProteinGroupLookup(CreateBaselineGroups());
        var groups = scores.Select(CreateTransientGroup).ToList();
        BaselineFdr.AlignProteinGroups(groupLookup, groups);

        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var psmLookup = BaselineFdr.CreatePsmLookup(CreateBaselinePsms(commonParameters));
        var psms = CreateTransientPsms(commonParameters, scores);
        BaselineFdr.AlignPsms(psmLookup, psms);

        Assert.Multiple(() =>
        {
            Assert.That(groups.Select(g => g.QValue), Is.EqualTo(expectedQValues).Within(1e-10));
            Assert.That(psms.Select(p => p.PsmFdrInfo.QValue), Is.EqualTo(expectedQValues).Within(1e-10));
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

    private static List<ProteinGroup> CreateBaselineGroups()
    {
        return
        [
            CreateProteinGroupWithFdr(score: 100, qValue: 0.01, bestPeptideQValue: 0.02, bestPeptidePep: 0.03, cumulativeTarget: 100, cumulativeDecoy: 1),
            CreateProteinGroupWithFdr(score: 50, qValue: 0.05, bestPeptideQValue: 0.06, bestPeptidePep: 0.07, cumulativeTarget: 80, cumulativeDecoy: 2),
            CreateProteinGroupWithFdr(score: 10, qValue: 0.10, bestPeptideQValue: 0.11, bestPeptidePep: 0.12, cumulativeTarget: 60, cumulativeDecoy: 3),
        ];
    }

    private static ProteinGroup CreateTransientGroup(double score)
    {
        return CreateProteinGroupWithFdr(score, qValue: 0.9, bestPeptideQValue: 0.9, bestPeptidePep: 0.9, cumulativeTarget: 1, cumulativeDecoy: 1);
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

    private static ProteinGroup CreateProteinGroupWithFdr(double score, double qValue, double bestPeptideQValue,
        double bestPeptidePep, int cumulativeTarget, int cumulativeDecoy)
    {
        var proteinGroup = ParallelSearchTestContextFactory.CreateProteinGroup(isDecoy: false, qValue: qValue, peptideCount: 1);
        proteinGroup.ProteinGroupScore = score;
        proteinGroup.BestPeptideQValue = bestPeptideQValue;
        proteinGroup.BestPeptidePEP = bestPeptidePep;
        proteinGroup.CumulativeTarget = cumulativeTarget;
        proteinGroup.CumulativeDecoy = cumulativeDecoy;
        return proteinGroup;
    }
}
