using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer.DatabaseLoading;
using TaskLayer;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;
using Readers;

namespace Test.Transcriptomics;

[TestFixture]
public class TestCalibration
{
    public static RnaSearchParameters SearchParameters;
    public static CommonParameters CommonParameters;

    [OneTimeSetUp]
    public static void Setup()
    {
        SearchParameters = new RnaSearchParameters
        {
            DecoyType = DecoyType.Reverse,
            MassDiffAcceptorType = MassDiffAcceptorType.Exact,
            DisposeOfFileWhenDone = true
        };
        CommonParameters = new CommonParameters
        (
            dissociationType: DissociationType.CID,
            deconvolutionMaxAssumedChargeState: -20,
            deconvolutionIntensityRatio: 3,
            deconvolutionMassTolerance: new PpmTolerance(20),
            precursorMassTolerance: new PpmTolerance(10),
            productMassTolerance: new PpmTolerance(5),
            scoreCutoff: 5,
            totalPartitions: 1,
            maxThreadsToUsePerFile: 1,
            doPrecursorDeconvolution: true,
            useProvidedPrecursorInfo: false,
            digestionParams: new RnaDigestionParams()
        );
    }

    [Test]
    public void CalibrationTask_ReducesMassError()
    {
        // Arrange
        string testDir = TestContext.CurrentContext.TestDirectory;
        List<string> dbPaths = 
            [Path.Combine(testDir, "Transcriptomics", "TestData", "16mer2.fasta")];
        List<string> spectraPaths = 
            [Path.Combine(testDir, "Transcriptomics", "TestData", "RnaStandard_Subset.mzML")];
        string outputDir = Path.Combine(testDir, "CalibrationTest");

        if (Directory.Exists(outputDir))
            Directory.Delete(outputDir, true); // Clean up previous test run
        Directory.CreateDirectory(outputDir);

        // Set up search-cali-search task for RNA
        var dbs = dbPaths.Select(dbPath => new DbForTask(dbPath, false)).ToList();
        var searchTask = new SearchTask()
        {
            SearchParameters = SearchParameters,
            CommonParameters = CommonParameters
        };
        var cali = new CalibrationTask()
        {
            CommonParameters = CommonParameters
        };
        var searchTask2 = new SearchTask()
        {
            SearchParameters = SearchParameters,
            CommonParameters = CommonParameters
        };
        var runner = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)>()
        {
            ("Search1", searchTask),
            ("Calibration", cali),
            ("Search2", searchTask2)
        }, spectraPaths, dbs, outputDir);
        runner.Run();

        var firstSearchDir = Path.Combine(outputDir, "Search1");
        var firstOsmPath = Directory.GetFiles(firstSearchDir, "*OSMs.osmtsv", SearchOption.AllDirectories).First();
        var firstOsms = SpectrumMatchTsvReader.ReadTsv(firstOsmPath, out var warnings);
        Assert.That(warnings.Count, Is.EqualTo(0));

        var secondSearchDir = Path.Combine(outputDir, "Search2");
        var secondOsmPath = Directory.GetFiles(secondSearchDir, "*OSMs.osmtsv", SearchOption.AllDirectories).First();
        var secondOsms = SpectrumMatchTsvReader.ReadTsv(secondOsmPath, out warnings);
        Assert.That(warnings.Count, Is.EqualTo(0));

        // Assert that calibration increased search results
        Assert.That(secondOsms.Count, Is.GreaterThan(firstOsms.Count), "Expected more OligoSpectralMatches after calibration.");
        Assert.That(secondOsms.Count(p => p.QValue <= 0.01), Is.GreaterThan(firstOsms.Count(p => p.QValue <= 0.01)),
            "Expected more OligoSpectralMatches with QValue <= 0.01 after calibration.");

        // Assert that the mass error is reduced after calibration
        var firstErrors = firstOsms.Select(p => p.MassDiffDa)
            .Select(double.Parse)
            .ToList();
        var secondErrors = secondOsms.Select(p => p.MassDiffDa)
            .Select(double.Parse)
            .ToList();

        var firstAverage = firstErrors.Average();
        var secondAverage = secondErrors.Average();

        Assert.That(secondAverage, Is.LessThan(firstAverage),
            "Expected average mass error to be reduced after calibration.");

        Directory.Delete(outputDir, true); // Clean up after test
    }
}
