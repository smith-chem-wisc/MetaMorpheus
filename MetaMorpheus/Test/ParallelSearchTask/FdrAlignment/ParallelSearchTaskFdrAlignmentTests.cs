using System;
using System.Collections.Generic;
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
    public void PsmSpectralMatchFdrAlignmentService_AlignsAndClamps()
    {
        var service = new PsmSpectralMatchFdrAlignmentService();

        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var baselinePsms = new List<SpectralMatch>
        {
            CreateSpectralMatchWithFdr(commonParameters, 100, psmQValue: 0.01, peptideQValue: 0.01, scanIndex: 101),
            CreateSpectralMatchWithFdr(commonParameters, 50, psmQValue: 0.05, peptideQValue: 0.05, scanIndex: 102),
            CreateSpectralMatchWithFdr(commonParameters, 10, psmQValue: 0.10, peptideQValue: 0.10, scanIndex: 103),
        };
        service.BuildBaselineCache(baselinePsms);

        var transientPsms = new List<SpectralMatch>
        {
            ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, false, 110, 0.9, 0.9, 1),
            ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, false, 80, 0.9, 0.9, 2),
            ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, false, 40, 0.9, 0.9, 3),
            ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, false, 5, 0.9, 0.9, 4),
        };

        var result = service.ApplyBaseline(transientPsms);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedCount, Is.EqualTo(4));
            Assert.That(result.ClampedHighCount, Is.EqualTo(1));
            Assert.That(result.ClampedLowCount, Is.EqualTo(1));
            Assert.That(transientPsms[0].PsmFdrInfo.QValue, Is.EqualTo(0.01).Within(1e-10));
            Assert.That(transientPsms[1].PsmFdrInfo.QValue, Is.EqualTo(0.01).Within(1e-10));
            Assert.That(transientPsms[2].PsmFdrInfo.QValue, Is.EqualTo(0.05).Within(1e-10));
            Assert.That(transientPsms[3].PsmFdrInfo.QValue, Is.EqualTo(0.10).Within(1e-10));
        });
    }

    [Test]
    public void PeptideSpectralMatchFdrAlignmentService_SkipsMissingPeptideFdr()
    {
        var service = new PeptideSpectralMatchFdrAlignmentService();

        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var baselinePsms = new List<SpectralMatch>
        {
            CreateSpectralMatchWithFdr(commonParameters, 100, psmQValue: 0.01, peptideQValue: null, scanIndex: 111),
            CreateSpectralMatchWithFdr(commonParameters, 50, psmQValue: 0.05, peptideQValue: 0.07, scanIndex: 112),
            CreateSpectralMatchWithFdr(commonParameters, 10, psmQValue: 0.10, peptideQValue: 0.2, scanIndex: 113),
        };
        service.BuildBaselineCache(baselinePsms);

        var transientPeptides = new List<SpectralMatch>
        {
            ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, false, 100, 0.9, 0.9, 11),
            ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, false, 40, 0.9, 0.9, 12),
        };

        var result = service.ApplyBaseline(transientPeptides);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedCount, Is.EqualTo(2));
            Assert.That(result.ClampedHighCount, Is.EqualTo(1));
            Assert.That(result.ClampedLowCount, Is.EqualTo(0));
            Assert.That(transientPeptides[0].PeptideFdrInfo.QValue, Is.EqualTo(0.07).Within(1e-10));
            Assert.That(transientPeptides[1].PeptideFdrInfo.QValue, Is.EqualTo(0.07).Within(1e-10));
        });
    }

    [Test]
    public void ProteinGroupFdrAlignmentService_AlignsAndClamps()
    {
        var service = new ProteinGroupFdrAlignmentService();

        var baselineGroups = new List<ProteinGroup>
        {
            CreateProteinGroupWithFdr(score: 100, qValue: 0.01, bestPeptideQValue: 0.02, bestPeptidePep: 0.03, cumulativeTarget: 100, cumulativeDecoy: 1),
            CreateProteinGroupWithFdr(score: 50, qValue: 0.05, bestPeptideQValue: 0.06, bestPeptidePep: 0.07, cumulativeTarget: 80, cumulativeDecoy: 2),
            CreateProteinGroupWithFdr(score: 10, qValue: 0.10, bestPeptideQValue: 0.11, bestPeptidePep: 0.12, cumulativeTarget: 60, cumulativeDecoy: 3),
        };
        service.BuildBaselineCache(baselineGroups);

        var transientGroups = new List<ProteinGroup>
        {
            CreateProteinGroupWithFdr(score: 110, qValue: 0.9, bestPeptideQValue: 0.9, bestPeptidePep: 0.9, cumulativeTarget: 1, cumulativeDecoy: 1),
            CreateProteinGroupWithFdr(score: 80, qValue: 0.9, bestPeptideQValue: 0.9, bestPeptidePep: 0.9, cumulativeTarget: 1, cumulativeDecoy: 1),
            CreateProteinGroupWithFdr(score: 40, qValue: 0.9, bestPeptideQValue: 0.9, bestPeptidePep: 0.9, cumulativeTarget: 1, cumulativeDecoy: 1),
            CreateProteinGroupWithFdr(score: 5, qValue: 0.9, bestPeptideQValue: 0.9, bestPeptidePep: 0.9, cumulativeTarget: 1, cumulativeDecoy: 1),
        };

        var result = service.ApplyBaseline(transientGroups);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedCount, Is.EqualTo(4));
            Assert.That(result.ClampedHighCount, Is.EqualTo(1));
            Assert.That(result.ClampedLowCount, Is.EqualTo(1));
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
