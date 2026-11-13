using Chromatography.RetentionTimeCalibration;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.Chronologer;
using EngineLayer;
using MzLibUtil;
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
    // The test factory's scans all come from this file.
    private const string TestFile = "TestFile";

    /// <summary>One minute per % acetonitrile, so an error in minutes equals observed minus Chronologer's raw output.</summary>
    private static readonly Dictionary<string, IrtCalibrationModel> Identity = new()
    {
        [TestFile] = IrtCalibration.Line((new RtMinutes(0), new Irt(0)), (new RtMinutes(1), new Irt(1)))
    };

    private static double PredictOneAtATime(SpectralMatch psm)
    {
        using var predictor = new ChronologerRetentionTimePredictor();
        var peptide = (IRetentionPredictable)psm.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer;
        return predictor.PredictRetentionTimeEquivalent(peptide, out _)!.Value;
    }

    [Test]
    public void CanCollectData_WithMissingLists_ReturnsFalse()
    {
        var context = new TransientDatabaseContext
        {
            AllPsms = null!,
            TransientPsms = new List<SpectralMatch>(),
            AllPeptides = new List<SpectralMatch>(),
            TransientPeptides = new List<SpectralMatch>(),
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
            transientPeptides: [highQValuePsm],
            retentionTimeCalibrations: Identity);

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
            transientPeptides: [psmLevelOnly],
            retentionTimeCalibrations: Identity);

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
            transientPeptides: [peptideLevelOnly],
            retentionTimeCalibrations: Identity);

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
            allPsms: [psm], transientPsms: [psm], allPeptides: [psm], transientPeptides: [psm], retentionTimeCalibrations: Identity);

        var result = new RetentionTimeCollector().CollectData(context);

        // The batched and single-peptide paths run different libtorch kernels, so they agree to float32 precision only.
        double expectedError = 40 - PredictOneAtATime(psm);
        Assert.Multiple(() =>
        {
            Assert.That((double[])result[RetentionTimeCollector.PsmAllRtErrors], Has.Length.EqualTo(1));
            Assert.That(((double[])result[RetentionTimeCollector.PsmAllRtErrors])[0], Is.EqualTo(expectedError).Within(1e-4));
            Assert.That(((double[])result[RetentionTimeCollector.PeptideAllRtErrors])[0], Is.EqualTo(expectedError).Within(1e-4));
        });
    }

    [Test]
    public void CollectData_ConvertsPredictionToMinutesThroughTheFilesCalibration()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var psm = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 100,
            psmQValue: 0.001, peptideQValue: 0.001, baseSequence: "LGEYGFQNALIVR", retentionTime: 40);

        // 2 % acetonitrile per minute, starting at 5 % at minute 0: minutes = (% ACN - 5) / 2
        var calibrations = new Dictionary<string, IrtCalibrationModel>
        {
            [TestFile] = IrtCalibration.Line((new RtMinutes(0), new Irt(5)), (new RtMinutes(10), new Irt(25)))
        };
        var context = ParallelSearchTestContextFactory.CreateContext(commonParameters,
            allPsms: [psm], transientPsms: [psm], allPeptides: [psm], transientPeptides: [psm], retentionTimeCalibrations: calibrations);

        var result = new RetentionTimeCollector().CollectData(context);

        double expectedError = 40 - (PredictOneAtATime(psm) - 5) / 2;
        Assert.Multiple(() =>
        {
            Assert.That(((double[])result[RetentionTimeCollector.PsmAllRtErrors])[0], Is.EqualTo(expectedError).Within(1e-4));
            Assert.That((double)result[RetentionTimeCollector.PsmMeanAbsoluteRtError], Is.EqualTo(Math.Abs(expectedError)).Within(1e-4));
        });
    }

    [Test]
    public void CollectData_FileWithoutCalibration_GetsNoErrors()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var psm = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 100,
            psmQValue: 0.001, peptideQValue: 0.001, baseSequence: "LGEYGFQNALIVR", retentionTime: 40);

        foreach (var calibrations in new[] { null, new Dictionary<string, IrtCalibrationModel>() })
        {
            var context = ParallelSearchTestContextFactory.CreateContext(commonParameters,
                allPsms: [psm], transientPsms: [psm], allPeptides: [psm], transientPeptides: [psm], retentionTimeCalibrations: calibrations);

            var result = new RetentionTimeCollector().CollectData(context);

            Assert.Multiple(() =>
            {
                Assert.That((double[])result[RetentionTimeCollector.PsmAllRtErrors], Is.Empty);
                Assert.That((double[])result[RetentionTimeCollector.PeptideAllRtErrors], Is.Empty);
                Assert.That(double.IsNaN((double)result[RetentionTimeCollector.PsmMeanAbsoluteRtError]), Is.True);
            });
        }
    }

    [Test]
    public void CollectData_SequenceTheModelCannotPredict_IsLeftOutWithoutFailing()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        // Chronologer refuses sequences shorter than seven residues.
        var psm = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 100,
            psmQValue: 0.001, peptideQValue: 0.001, baseSequence: "PEPK", retentionTime: 10);
        var context = ParallelSearchTestContextFactory.CreateContext(commonParameters,
            allPsms: [psm], transientPsms: [psm], allPeptides: [psm], transientPeptides: [psm], retentionTimeCalibrations: Identity);

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
                allPsms: psms, transientPsms: psms, allPeptides: psms, transientPeptides: psms, retentionTimeCalibrations: Identity);
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

    /// <summary>Distinct tryptic-looking sequences (no internal K or R), long enough for Chronologer.</summary>
    private static List<string> SyntheticPeptides(int count, int seed)
    {
        const string residues = "ACDEFGHILMNQSTVWY";
        var random = new Random(seed);
        var peptides = new HashSet<string>();
        while (peptides.Count < count)
        {
            int length = random.Next(8, 20);
            peptides.Add(new string(Enumerable.Range(0, length).Select(_ => residues[random.Next(residues.Length)]).ToArray()) + "K");
        }
        return peptides.OrderBy(p => p, StringComparer.Ordinal).ToList();
    }

    [Test]
    public void FitRetentionTimeCalibrations_RecoversTheLineFromPredictionToMinutes()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        using var predictor = new ChronologerRetentionTimePredictor();

        // Observed minutes laid exactly on minutes = 0.5 * %ACN + 3, so the fit has nothing to smooth.
        var baseline = new List<SpectralMatch>();
        foreach (var (sequence, i) in SyntheticPeptides(40, seed: 7).Select((s, i) => (s, i)))
        {
            var probe = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, false, 100, 0.001, 0.001, baseSequence: sequence);
            double minutes = 0.5 * PredictOneAtATime(probe) + 3;
            baseline.Add(ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 100,
                psmQValue: 0.001, peptideQValue: 0.001, scanNumber: i + 1, baseSequence: sequence, retentionTime: minutes));
        }
        // Not confident, decoy: never anchors, however far off their retention times are.
        baseline.Add(ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, false, 100, 0.5, 0.5, scanNumber: 100,
            baseSequence: "LGEYGFQNALIVR", retentionTime: 500));
        baseline.Add(ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, true, 100, 0.001, 0.001, scanNumber: 101,
            baseSequence: "LGEYGFQNALIVR", retentionTime: 500));

        var warnings = new List<string>();
        var calibrations = RetentionTimeCollector.FitRetentionTimeCalibrations(baseline, 0.01, 2, warnings.Add);

        Assert.That(warnings, Is.Empty);
        Assert.That(calibrations.Keys, Is.EquivalentTo(new[] { TestFile }));
        var calibration = calibrations[TestFile];
        Assert.Multiple(() =>
        {
            // At most the 40 confident targets; the fit may set aside a few in sparse stretches.
            Assert.That(calibration.AnchorCount, Is.InRange(20, 40));
            foreach (double acetonitrile in new[] { 10.0, 20.0, 30.0 })
                Assert.That(calibration.ToRtMinutes(new Irt(acetonitrile)).Value, Is.EqualTo(0.5 * acetonitrile + 3).Within(0.01));
        });
    }

    [Test]
    public void FitRetentionTimeCalibrations_WithTooFewAnchors_WarnsAndLeavesTheFileOut()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var baseline = SyntheticPeptides(5, seed: 11)
            .Select((sequence, i) => ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 100,
                psmQValue: 0.001, peptideQValue: 0.001, scanNumber: i + 1, baseSequence: sequence, retentionTime: 10 + i))
            .ToList();

        var warnings = new List<string>();
        var calibrations = RetentionTimeCollector.FitRetentionTimeCalibrations(baseline, 0.01, 1, warnings.Add);

        Assert.Multiple(() =>
        {
            Assert.That(calibrations, Is.Empty);
            Assert.That(warnings, Has.Count.EqualTo(1));
            Assert.That(warnings[0], Does.Contain(TestFile));
        });
    }
}
