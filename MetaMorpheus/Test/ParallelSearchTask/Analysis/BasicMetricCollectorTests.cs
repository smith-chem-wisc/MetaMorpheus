using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using TaskLayer.ParallelSearch.Analysis;
using TaskLayer.ParallelSearch.Analysis.Collectors;
using Test.ParallelSearchTask.Utility;

namespace Test.ParallelSearchTask.Analysis;

[TestFixture]
public class BasicMetricCollectorTests
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
            TransientDatabase = null!,
            TransientProteins = null!,
            TransientProteinAccessions = new HashSet<string>(),
        };

        var collector = new BasicMetricCollector();

        Assert.That(collector.CanCollectData(context), Is.False);
    }

    [Test]
    public void CollectData_UsesMinimumThresholdAcrossPsmAndPeptide()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters(qValueThreshold: 0.01, pepQValueThreshold: 0.05);

        var allPsms = new List<EngineLayer.SpectralMatch>
        {
            ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 100, psmQValue: 0.009, peptideQValue: 0.009),
            ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 90, psmQValue: 0.02, peptideQValue: 0.02),
            ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: true, score: 80, psmQValue: 0.001, peptideQValue: 0.001),
        };

        var allPeptides = new List<EngineLayer.SpectralMatch>
        {
            ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 75, psmQValue: 0.5, peptideQValue: 0.009),
            ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 65, psmQValue: 0.5, peptideQValue: 0.02),
            ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: true, score: 55, psmQValue: 0.5, peptideQValue: 0.001),
        };

        var context = ParallelSearchTestContextFactory.CreateContext(
            commonParameters,
            allPsms,
            transientPsms: [allPsms[0], allPsms[2]],
            allPeptides,
            transientPeptides: [allPeptides[0], allPeptides[2]],
            totalProteins: 10,
            transientPeptideCount: 4);

        var collector = new BasicMetricCollector();
        var result = collector.CollectData(context);

        Assert.Multiple(() =>
        {
            Assert.That(result[BasicMetricCollector.TotalProteins], Is.EqualTo(10));
            Assert.That(result[BasicMetricCollector.TransientProteinCount], Is.EqualTo(2));
            Assert.That(result[BasicMetricCollector.TransientPeptideCount], Is.EqualTo(4));
            Assert.That(result[BasicMetricCollector.TargetPsmsAtQValueThreshold], Is.EqualTo(1));
            Assert.That(result[BasicMetricCollector.TargetPsmsFromTransientDb], Is.EqualTo(1));
            Assert.That(result[BasicMetricCollector.TargetPsmsFromTransientDbAtQValueThreshold], Is.EqualTo(1));
            Assert.That(result[BasicMetricCollector.TargetPeptidesAtQValueThreshold], Is.EqualTo(1));
            Assert.That(result[BasicMetricCollector.TargetPeptidesFromTransientDb], Is.EqualTo(1));
            Assert.That(result[BasicMetricCollector.TargetPeptidesFromTransientDbAtQValueThreshold], Is.EqualTo(1));
        });
    }

    [Test]
    [TestCase(120, 60)] // PEP filtering: counted by PEP q-value at the PEP q-value threshold, as FilteredPsms does
    [TestCase(50, 50)]  // fewer than 100 matches: FilteredPsms does not filter, so every match counts
    public void CollectData_WithAPepQValueThresholdBelowTheQValueThreshold_FiltersAsFilteredPsmsDoes(int matchCount, int expectedConfident)
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters(qValueThreshold: 1.0, pepQValueThreshold: 0.01);
        // q-values all above 0.01, PEP q-values passing for the first half only
        var psms = Enumerable.Range(0, matchCount).Select(i =>
        {
            var psm = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 10,
                psmQValue: 0.5, peptideQValue: 0.5, scanNumber: i + 1);
            double pepQValue = i < matchCount / 2 ? 0.001 : 0.5;
            psm.SetFdrValues(0, 0, 0.5, 0, 0, 0.5, 0, pepQValue);
            psm.PeptideFdrInfo.PEP_QValue = pepQValue;
            return psm;
        }).ToList();
        var context = ParallelSearchTestContextFactory.CreateContext(commonParameters, allPsms: psms, transientPsms: psms,
            allPeptides: psms, transientPeptides: psms);

        var result = new BasicMetricCollector().CollectData(context);

        Assert.Multiple(() =>
        {
            Assert.That(result[BasicMetricCollector.TargetPsmsAtQValueThreshold], Is.EqualTo(expectedConfident));
            Assert.That(result[BasicMetricCollector.TargetPsmsFromTransientDbAtQValueThreshold], Is.EqualTo(expectedConfident));
            Assert.That(result[BasicMetricCollector.TargetPeptidesFromTransientDbAtQValueThreshold], Is.EqualTo(expectedConfident));
        });
    }

    [Test]
    public void CollectData_OnTheQValuePath_AlsoHoldsTheNotchQValueToTheThreshold()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        // every q-value passes, but the notch q-value passes for the first half only
        var psms = Enumerable.Range(0, 150).Select(i =>
        {
            var psm = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 10,
                psmQValue: 0.005, peptideQValue: 0.005, scanNumber: i + 1);
            double notch = i < 75 ? 0.005 : 0.05;
            psm.SetFdrValues(0, 0, 0.005, 0, 0, notch, 0, 0.005);
            psm.PeptideFdrInfo.QValueNotch = notch;
            return psm;
        }).ToList();
        var context = ParallelSearchTestContextFactory.CreateContext(commonParameters, allPsms: psms, transientPsms: psms,
            allPeptides: psms, transientPeptides: psms);

        var result = new BasicMetricCollector().CollectData(context);

        Assert.Multiple(() =>
        {
            Assert.That(result[BasicMetricCollector.TargetPsmsAtQValueThreshold], Is.EqualTo(75));
            Assert.That(result[BasicMetricCollector.TargetPsmsFromTransientDbAtQValueThreshold], Is.EqualTo(75));
            Assert.That(result[BasicMetricCollector.TargetPeptidesFromTransientDbAtQValueThreshold], Is.EqualTo(75));
        });
    }
}
