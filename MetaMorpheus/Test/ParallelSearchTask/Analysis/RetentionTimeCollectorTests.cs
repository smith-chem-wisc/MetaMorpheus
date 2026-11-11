using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.Chronologer;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
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

    [Test]
    public void CollectData_BatchedPredictionMatchesSinglePeptidePrediction()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var psm = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 100,
            psmQValue: 0.001, peptideQValue: 0.001, baseSequence: "LGEYGFQNALIVR", retentionTime: 40);
        var context = ParallelSearchTestContextFactory.CreateContext(commonParameters,
            allPsms: [psm], transientPsms: [psm], allPeptides: [psm], transientPeptides: [psm]);

        var result = new RetentionTimeCollector().CollectData(context);

        // The batched and single-peptide paths run different libtorch kernels, so they agree to float32 precision only.
        using var predictor = new ChronologerRetentionTimePredictor();
        var peptide = (IRetentionPredictable)psm.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer;
        double expectedError = 40 - predictor.PredictRetentionTimeEquivalent(peptide, out _)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That((double[])result[RetentionTimeCollector.PsmAllRtErrors], Has.Length.EqualTo(1));
            Assert.That(((double[])result[RetentionTimeCollector.PsmAllRtErrors])[0], Is.EqualTo(expectedError).Within(1e-4));
            Assert.That(((double[])result[RetentionTimeCollector.PeptideAllRtErrors])[0], Is.EqualTo(expectedError).Within(1e-4));
        });
    }

    [Test]
    public void CollectData_SequenceTheModelCannotPredict_IsLeftOutWithoutFailing()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        // Chronologer refuses sequences shorter than seven residues.
        var psm = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 100,
            psmQValue: 0.001, peptideQValue: 0.001, baseSequence: "PEPK", retentionTime: 10);
        var context = ParallelSearchTestContextFactory.CreateContext(commonParameters,
            allPsms: [psm], transientPsms: [psm], allPeptides: [psm], transientPeptides: [psm]);

        var collector = new RetentionTimeCollector();
        var first = collector.CollectData(context);
        var second = collector.CollectData(context); // now answered from the cached failure

        Assert.Multiple(() =>
        {
            foreach (var result in new[] { first, second })
            {
                Assert.That((double[])result[RetentionTimeCollector.PsmAllRtErrors], Is.Empty);
                Assert.That(double.IsNaN((double)result[RetentionTimeCollector.PsmMeanAbsoluteRtError]), Is.True);
            }
        });
    }

    [Test]
    public void CollectData_DatabasesCollectedAtOnce_GetTheSameValueForASharedSequence()
    {
        // Sequences no other test predicts, so the shared cache starts without them.
        string[] sequences = ["AEFVEVTK", "ADLAGHGVSPLIK", "TVAAPSVFIFPPSDEQLK", "VLDALDSIK", "GDLGIEIPAEK", "LVNELTEFAK"];
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();

        // Four databases, each holding an overlapping subset in a different order, as transient databases would.
        var databases = Enumerable.Range(0, 4)
            .Select(d => sequences.Skip(d).Concat(sequences.Take(d)).Take(4)
                .Select((sequence, i) => ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false,
                    score: 100, psmQValue: 0.001, peptideQValue: 0.001, scanNumber: i + 1, baseSequence: sequence, retentionTime: 30))
                .ToList())
            .ToList();

        var errorsBySequence = new Dictionary<string, double>[databases.Count];
        Parallel.For(0, databases.Count, new ParallelOptions { MaxDegreeOfParallelism = databases.Count }, d =>
        {
            var psms = databases[d];
            var context = ParallelSearchTestContextFactory.CreateContext(commonParameters,
                allPsms: psms, transientPsms: psms, allPeptides: psms, transientPeptides: psms);
            var errors = (double[])new RetentionTimeCollector().CollectData(context)[RetentionTimeCollector.PsmAllRtErrors];
            errorsBySequence[d] = psms.Select(p => p.BaseSequence).Zip(errors).ToDictionary(t => t.First, t => t.Second);
        });

        // Bitwise equal: a batch's own value may differ from another batch's in the last float32 digits,
        // so every database must read the one value that reached the cache first.
        foreach (string sequence in sequences)
        {
            var values = errorsBySequence.Where(e => e.ContainsKey(sequence)).Select(e => e[sequence]).Distinct().ToList();
            Assert.That(values, Has.Count.EqualTo(1), sequence);
        }
    }
}
