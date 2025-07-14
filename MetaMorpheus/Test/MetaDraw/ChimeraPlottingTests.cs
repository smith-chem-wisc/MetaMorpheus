using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EngineLayer;
using EngineLayer.FdrAnalysis;
using GuiFunctions;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using OxyPlot;
using Readers;
using Path = System.IO.Path;

namespace Test.MetaDraw;

[TestFixture]
[ExcludeFromCodeCoverage]
public class ChimeraPlottingTests
{
    public static ChimeraTestCase OneProteinTwoProteoformChimeraGroup;
    public static ChimeraTestCase TwoProteinsTwoProteoformChimeraGroup;
    public static List<SpectrumMatchFromTsv> AllMatches;
    public static MsDataFile DataFile;
    public static string TestExportDirectory => Path.Combine(TestContext.CurrentContext.TestDirectory, "MetaDraw", "ChimeraPlottingTests");

    public record ChimeraTestCase(ChimeraGroupViewModel ChimeraGroup, Dictionary<OxyColor, List<MatchedFragmentIon>> ExpectedIonsByColor);

    static ChimeraPlottingTests()
    {
        var psmPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TDGPTMDSearchResults.psmtsv");
        string msDataPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TDGPTMDSearchSingleSpectra.mzML");

        DataFile = MsDataFileReader.GetDataFile(msDataPath).LoadAllStaticData();
        DataFile.InitiateDynamicConnection();
        AllMatches = SpectrumMatchTsvReader.ReadTsv(psmPath, out var warnings).Where(p => p.Ms2ScanNumber <= DataFile.Scans.Length).ToList();

        var testMs1Scan = DataFile.GetOneBasedScan(900);
        var testMs2Scan = DataFile.GetOneBasedScan(901);
        var testPsms = AllMatches.Where(p => p.Ms2ScanNumber == testMs2Scan.OneBasedScanNumber);
        var group = new ChimeraGroupViewModel(testPsms, testMs1Scan, testMs2Scan);
        var ions = group.AssignFragmentIonColors()
            .ToDictionary(p => p.Key, p => p.Value.Select(m => m.Item1).ToList());
        OneProteinTwoProteoformChimeraGroup = new ChimeraTestCase(group, ions);

        var testMs1Scan2 = DataFile.GetOneBasedScan(1243);
        var testMs2Scan2 = DataFile.GetOneBasedScan(1246);
        var testPsms2 = AllMatches.Where(p => p.Ms2ScanNumber == testMs2Scan2.OneBasedScanNumber && p.QValue <= 0.01);
        var group2 = new ChimeraGroupViewModel(testPsms2, testMs1Scan2, testMs2Scan2);
        var ions2 = group2.AssignFragmentIonColors()
            .ToDictionary(p => p.Key, p => p.Value.Select(m => m.Item1).ToList());

        TwoProteinsTwoProteoformChimeraGroup = new ChimeraTestCase(group2, ions2);
    }

    [OneTimeSetUp]
    public static void OneTimeSetup()
    {
        GlobalVariables.AnalyteType = AnalyteType.Proteoform;
        // Ensure the export directory exists in a new state
        if (Directory.Exists(TestExportDirectory))
        {
            Directory.Delete(TestExportDirectory, true);
        }
        Directory.CreateDirectory(TestExportDirectory);
    }

    [OneTimeTearDown]
    public static void OneTimeTearDown()
    {
        GlobalVariables.AnalyteType = AnalyteType.Peptide;
        DataFile.CloseDynamicConnection();
        // Clean up the export directory after tests
        if (Directory.Exists(TestExportDirectory))
        {
            Directory.Delete(TestExportDirectory, true);
        }
    }

    #region Tab View Model

    [Test]
    [NonParallelizable]
    public void ChimeraTabViewModel_ResultsFilteredCorrectly()
    {
        var dataFiles = new Dictionary<string, MsDataFile>()
        {
            {"FXN3_tr1_032017-calib", DataFile }
        };

        var chimeraAnalysisTabViewModel = new ChimeraAnalysisTabViewModel(AllMatches, dataFiles, TestExportDirectory);
        Assert.That(chimeraAnalysisTabViewModel.ChimeraGroupViewModels.All(p => p.Count > 1),
            Is.True, "All chimera groups should have at least two PSMs.");
        Assert.That(chimeraAnalysisTabViewModel.ChimeraGroupViewModels.Count, Is.GreaterThan(0), "There should be at least one chimera group.");

        // Check that all PSMs in chimera groups pass the Q-value filter and do not contain decoys if ShowDecoys is false
        chimeraAnalysisTabViewModel.ChimeraGroupViewModels.SelectMany(p => p.ChimericPsms)
            .ToList()
            .ForEach(p =>
            {
                Assert.That(p.Psm.QValue, Is.LessThanOrEqualTo(MetaDrawSettings.QValueFilter), "Chimeric PSMs should pass the Q-value filter.");
                Assert.That(p.Psm.DecoyContamTarget, Does.Not.Contain('D'), "Chimeric PSMs should not contain decoys if ShowDecoys is false.");
            });

        MetaDrawSettings.QValueFilter = 0.5; // set filter higher
        MetaDrawSettings.ShowDecoys = true; // show decoys

        chimeraAnalysisTabViewModel = new ChimeraAnalysisTabViewModel(AllMatches, dataFiles, TestExportDirectory);
        Assert.That(chimeraAnalysisTabViewModel.ChimeraGroupViewModels.All(p => p.Count > 1),
            Is.True, "All chimera groups should have at least two PSMs after changing settings.");
        Assert.That(chimeraAnalysisTabViewModel.ChimeraGroupViewModels.Count, Is.GreaterThan(0), "There should still be at least one chimera group after changing settings.");

        // Check that all PSMs in chimera groups pass the new Q-value filter and may contain decoys
        chimeraAnalysisTabViewModel.ChimeraGroupViewModels.SelectMany(p => p.ChimericPsms)
            .ToList()
            .ForEach(p =>
            {
                Assert.That(p.Psm.QValue, Is.LessThanOrEqualTo(MetaDrawSettings.QValueFilter), "Chimeric PSMs should pass the new Q-value filter.");
            });
        Assert.That(chimeraAnalysisTabViewModel.ChimeraGroupViewModels.Any(p => p.ChimericPsms.Any(q => q.Psm.DecoyContamTarget.Contains('D'))),
            Is.True, "Some chimeric PSMs should contain decoys after changing ShowDecoys to true.");

        // Reset settings for other tests
        MetaDrawSettings.QValueFilter = 0.01; // reset filter
        MetaDrawSettings.ShowDecoys = false; // reset to not show decoys
    }

    [Test]
    public static void ChimeraTabViewModel_PrecursorAssignmentIsCorrectInGroups()
    {
        var dataFiles = new Dictionary<string, MsDataFile>()
        {
            {"FXN3_tr1_032017-calib", DataFile }
        };
        var chimeraAnalysisTabViewModel = new ChimeraAnalysisTabViewModel(AllMatches, dataFiles, TestExportDirectory);

        foreach (var chimeraGroup in chimeraAnalysisTabViewModel.ChimeraGroupViewModels)
        {
            // Only test groups with at least 2 PSMs
            if (chimeraGroup.Count < 2)
                continue;

            var chimericPsms = chimeraGroup.ChimericPsms.ToList();
            for (int i = 0; i < chimericPsms.Count; i++)
            {
                var psm = chimericPsms[i].Psm;
                var envelope = chimericPsms[i].PrecursorEnvelope;

                // Compare to all other envelopes in the group
                double minDiff = Math.Abs(psm.PrecursorMass - envelope.MonoisotopicMass);
                for (int j = 0; j < chimericPsms.Count; j++)
                {
                    if (i == j) continue;
                    var otherEnvelope = chimericPsms[j].PrecursorEnvelope;
                    double diff = Math.Abs(psm.PrecursorMass - otherEnvelope.MonoisotopicMass);
                    Assert.That(minDiff, Is.LessThanOrEqualTo(diff),
                        $"PSM {i} in group (scan {chimeraGroup.Ms2ScanNumber}) should be assigned to its closest envelope.");
                }
            }
        }
    }

    #endregion

    #region Group View Model

    [Test]
    public static void ChimeraGroupViewModel_FragmentCounts_OneProteinTwoProteoforms()
    {
        // Arrange
        var chimeraGroup = OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        Assert.That(chimeraGroup.Count, Is.EqualTo(2));
        Assert.That(chimeraGroup.ProteinCount, Is.EqualTo(1));

        var terminalFragments = chimeraGroup.ChimericPsms.SelectMany(p => p.Psm.MatchedIons)
            .Where(p => !p.IsInternalFragment)
            .ToList();
        Assert.That(terminalFragments.Count, Is.EqualTo(chimeraGroup.TotalFragments),  chimeraGroup.TotalFragments + " terminal fragments should match total fragments in the group.");

        var uniqueFragments = terminalFragments
            .Distinct()
            .ToList();
        Assert.That(uniqueFragments.Count, Is.EqualTo(chimeraGroup.UniqueFragments), "Unique fragments should match total unique fragments in the group.");
    }

    [Test]
    public static void ChimeraGroupViewModel_FragmentCounts_TwoProteinsTwoProteoforms()
    {
        // Arrange
        var chimeraGroup = TwoProteinsTwoProteoformChimeraGroup.ChimeraGroup;
        Assert.That(chimeraGroup.Count, Is.EqualTo(2));
        Assert.That(chimeraGroup.ProteinCount, Is.EqualTo(2));

        var terminalFragments = chimeraGroup.ChimericPsms.SelectMany(p => p.Psm.MatchedIons)
            .Where(p => !p.IsInternalFragment)
            .ToList();
        Assert.That(terminalFragments.Count, Is.EqualTo(chimeraGroup.TotalFragments), chimeraGroup.TotalFragments + " terminal fragments should match total fragments in the group.");

        var uniqueFragments = terminalFragments
            .Distinct()
            .ToList();
        Assert.That(uniqueFragments.Count, Is.EqualTo(chimeraGroup.UniqueFragments), "Unique fragments should match total unique fragments in the group.");
    }

    [Test]
    public static void ChimeraGroupViewModel_PrecursorAssignmentIsCorrect_OneProteinTwoProteoforms()
    {
        // Arrange
        var chimeraGroup = OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        Assert.That(chimeraGroup.Count, Is.EqualTo(2));
        Assert.That(chimeraGroup.ProteinCount, Is.EqualTo(1));


        var envelope1 = chimeraGroup.ChimericPsms[0].PrecursorEnvelope;
        var envelope2 = chimeraGroup.ChimericPsms[1].PrecursorEnvelope;
        var psm1 = chimeraGroup.ChimericPsms[0].Psm;
        var psm2 = chimeraGroup.ChimericPsms[1].Psm;


        // ensure mass differences are minimized
        var psm1ToEnvelope1 = Math.Abs(psm1.PrecursorMass - envelope1.MonoisotopicMass);
        var psm2ToEnvelope2 = Math.Abs(psm2.PrecursorMass - envelope2.MonoisotopicMass);
        var psm1ToEnvelope2 = Math.Abs(psm1.PrecursorMass - envelope2.MonoisotopicMass);
        var psm2ToEnvelope1 = Math.Abs(psm2.PrecursorMass - envelope1.MonoisotopicMass);

        Assert.That(psm1ToEnvelope1, Is.LessThanOrEqualTo(psm1ToEnvelope2), "PSM 1 should be assigned to its own envelope.");
        Assert.That(psm2ToEnvelope2, Is.LessThanOrEqualTo(psm2ToEnvelope1), "PSM 2 should be assigned to its own envelope.");
    }

    [Test]
    public static void ChimeraGroupViewModel_PrecursorAssignmentIsCorrect_TwoProteinsTwoProteoforms()
    {
        // Arrange
        var chimeraGroup = TwoProteinsTwoProteoformChimeraGroup.ChimeraGroup;
        Assert.That(chimeraGroup.Count, Is.EqualTo(2));
        Assert.That(chimeraGroup.ProteinCount, Is.EqualTo(2));

        var envelope1 = chimeraGroup.ChimericPsms[0].PrecursorEnvelope;
        var envelope2 = chimeraGroup.ChimericPsms[1].PrecursorEnvelope;
        var psm1 = chimeraGroup.ChimericPsms[0].Psm;
        var psm2 = chimeraGroup.ChimericPsms[1].Psm;


        // ensure mass differences are minimized
        var psm1ToEnvelope1 = Math.Abs(psm1.PrecursorMass - envelope1.MonoisotopicMass);
        var psm2ToEnvelope2 = Math.Abs(psm2.PrecursorMass - envelope2.MonoisotopicMass);
        var psm1ToEnvelope2 = Math.Abs(psm1.PrecursorMass - envelope2.MonoisotopicMass);
        var psm2ToEnvelope1 = Math.Abs(psm2.PrecursorMass - envelope1.MonoisotopicMass);

        Assert.That(psm1ToEnvelope1, Is.LessThanOrEqualTo(psm1ToEnvelope2), "PSM 1 should be assigned to its own envelope.");
        Assert.That(psm2ToEnvelope2, Is.LessThanOrEqualTo(psm2ToEnvelope1), "PSM 2 should be assigned to its own envelope.");
    }

    [Test]
    public static void ChimeraGroupViewModel_PrecursorColorAssignmentIsCorrect_OneProteinTwoProteoforms()
    {
        // Arrange
        var chimeraGroup = OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        Assert.That(chimeraGroup.Count, Is.EqualTo(2));

        var envelope1 = chimeraGroup.ChimericPsms[0].PrecursorEnvelope;
        var envelope2 = chimeraGroup.ChimericPsms[1].PrecursorEnvelope;

        var envelope1Peaks = envelope1.Peaks;
        var envelope2Peaks = envelope2.Peaks;

        var sharedPeakColor = ChimeraSpectrumMatchPlot.ColorByProteinDictionary[0][0];
        var firstColor = ChimeraSpectrumMatchPlot.ColorByProteinDictionary[0][1];
        var secondColor = ChimeraSpectrumMatchPlot.ColorByProteinDictionary[0][2];

        // One color for each of the precursor envelopes. 
        Assert.That(chimeraGroup.PrecursorIonsByColor.Count(), Is.EqualTo(2));
        Assert.That(chimeraGroup.PrecursorIonsByColor.ContainsKey(sharedPeakColor), Is.False, "Shared peak color should not be present.");
        Assert.That(chimeraGroup.PrecursorIonsByColor.ContainsKey(firstColor), Is.True, "First color should be present.");
        Assert.That(chimeraGroup.PrecursorIonsByColor.ContainsKey(secondColor), Is.True, "Second color should be present.");

        Assert.That(chimeraGroup.PrecursorIonsByColor[firstColor].Count, Is.EqualTo(envelope1Peaks.Count), "First color should match envelope 1 peaks.");
        Assert.That(chimeraGroup.PrecursorIonsByColor[secondColor].Count, Is.EqualTo(envelope2Peaks.Count), "Second color should match envelope 2 peaks.");
        Assert.That(chimeraGroup.ChimericPsms.All(p => p.ProteinColor == sharedPeakColor));
    }

    [Test]
    public static void ChimeraGroupViewModel_PrecursorColorAssignmentIsCorrect_TwoProteinsTwoProteoforms()
    {
        // Arrange
        var chimeraGroup = TwoProteinsTwoProteoformChimeraGroup.ChimeraGroup;
        Assert.That(chimeraGroup.Count, Is.EqualTo(2));

        var envelope1 = chimeraGroup.ChimericPsms[0].PrecursorEnvelope;
        var envelope2 = chimeraGroup.ChimericPsms[1].PrecursorEnvelope;

        var envelope1Peaks = envelope1.Peaks;
        var envelope2Peaks = envelope2.Peaks;

        var sharedPeakColor = OxyColors.Black;
        var firstProteinColor = ChimeraSpectrumMatchPlot.ColorByProteinDictionary[0][0];
        var firstProteoformOfFirstProteinColor = ChimeraSpectrumMatchPlot.ColorByProteinDictionary[0][1];
        var secondProteinColor = ChimeraSpectrumMatchPlot.ColorByProteinDictionary[1][0];
        var firstProteoformOfSecondProteinColor = ChimeraSpectrumMatchPlot.ColorByProteinDictionary[1][1];

        // One color for each of the precursor envelopes. 
        Assert.That(chimeraGroup.PrecursorIonsByColor.Count(), Is.EqualTo(2));
        Assert.That(chimeraGroup.PrecursorIonsByColor.ContainsKey(sharedPeakColor), Is.False, "Shared peak color should not be present.");
        Assert.That(chimeraGroup.PrecursorIonsByColor.ContainsKey(firstProteoformOfFirstProteinColor), Is.True, "First color should be present.");
        Assert.That(chimeraGroup.PrecursorIonsByColor.ContainsKey(firstProteoformOfSecondProteinColor), Is.True, "Second color should be present.");

        Assert.That(chimeraGroup.PrecursorIonsByColor[firstProteoformOfFirstProteinColor].Count, Is.EqualTo(envelope1Peaks.Count), "First color should match envelope 1 peaks.");
        Assert.That(chimeraGroup.PrecursorIonsByColor[firstProteoformOfSecondProteinColor].Count, Is.EqualTo(envelope2Peaks.Count), "Second color should match envelope 2 peaks.");
        Assert.That(chimeraGroup.ChimericPsms.Any(p => p.ProteinColor == firstProteinColor));
        Assert.That(chimeraGroup.ChimericPsms.Any(p => p.Color == firstProteoformOfFirstProteinColor));
        Assert.That(chimeraGroup.ChimericPsms.Any(p => p.ProteinColor == secondProteinColor));
        Assert.That(chimeraGroup.ChimericPsms.Any(p => p.Color == firstProteoformOfSecondProteinColor));
    }

    [Test]
    public static void ChimeraGroupViewModel_FragmentColorAssignmentIsCorrect_OneProteinTwoProteoforms()
    {
        // Arrange
        var chimeraGroup = OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        var expectedIons = OneProteinTwoProteoformChimeraGroup.ExpectedIonsByColor;

        // Assert that the expected ions match the actual ions in the group
        Assert.That(chimeraGroup.MatchedFragmentIonsByColor.Count, Is.EqualTo(expectedIons.Count), "Precursor ions by color count should match expected ions count.");

        foreach (var ionsByColor in chimeraGroup.MatchedFragmentIonsByColor)
        {
            var color = ionsByColor.Key;
            var actualIons = ionsByColor.Value;
            Assert.That(expectedIons.ContainsKey(color), Is.True, $"Expected ions should contain color {color}.");

            var expectedIonsList = expectedIons[color];
            // Check that the actual ions match the expected ions
            Assert.That(actualIons.Count, Is.EqualTo(expectedIonsList.Count), $"Count of matched fragment ions for color {color} should match expected count.");
            foreach (var ion in actualIons)
            {
                Assert.That(expectedIonsList.Any(e => Math.Abs(e.Mz - ion.Item1.Mz) < 1e-3 && e.NeutralTheoreticalProduct.ProductType == ion.Item1.NeutralTheoreticalProduct.ProductType), Is.True,
                    $"Matched fragment ion {ion} should be in expected ions for color {color}.");
            }

        }
    }

    [Test]
    public static void ChimeraGroupViewModel_FragmentColorAssignmentIsCorrect_TwoProteinsTwoProteoforms()
    {
        // Arrange
        var chimeraGroup = TwoProteinsTwoProteoformChimeraGroup.ChimeraGroup;
        var expectedIons = TwoProteinsTwoProteoformChimeraGroup.ExpectedIonsByColor;

        // Assert that the expected ions match the actual ions in the group
        Assert.That(chimeraGroup.MatchedFragmentIonsByColor.Count, Is.EqualTo(expectedIons.Count), "Precursor ions by color count should match expected ions count.");

        foreach (var ionsByColor in chimeraGroup.MatchedFragmentIonsByColor)
        {
            var color = ionsByColor.Key;
            var actualIons = ionsByColor.Value;
            Assert.That(expectedIons.ContainsKey(color), Is.True, $"Expected ions should contain color {color}.");

            var expectedIonsList = expectedIons[color];
            // Check that the actual ions match the expected ions
            Assert.That(actualIons.Count, Is.EqualTo(expectedIonsList.Count), $"Count of matched fragment ions for color {color} should match expected count.");
            foreach (var ion in actualIons)
            {
                Assert.That(expectedIonsList.Any(e => Math.Abs(e.Mz - ion.Item1.Mz) < 1e-3 && e.NeutralTheoreticalProduct.ProductType == ion.Item1.NeutralTheoreticalProduct.ProductType), Is.True,
                    $"Matched fragment ion {ion} should be in expected ions for color {color}.");
            }

        }
    }


    #endregion
}

[TestFixture]
public class HungarianAlgorithmTests
{
    [Test]
    public void SimpleSquareMatrix_AssignsDiagonal()
    {
        // Arrange
        double[,] cost = new double[,]
        {
            { 1, 2, 3 },
            { 2, 1, 3 },
            { 3, 2, 1 }
        };

        // Act
        int[] assignment = HungarianAlgorithm.FindAssignments(cost);

        // Assert
        Assert.That(assignment, Is.EqualTo(new[] { 0, 1, 2 }));
    }

    [Test]
    public void RectangularMatrix_AssignsOptimal()
    {
        // Arrange: 2x3 matrix, optimal is row 0 to col 1, row 1 to col 2
        double[,] cost = new double[,]
        {
            { 10, 1, 10 },
            { 10, 10, 1 }
        };

        // Act
        int[] assignment = HungarianAlgorithm.FindAssignments(cost);

        // Assert
        Assert.That(assignment[0], Is.EqualTo(1));
        Assert.That(assignment[1], Is.EqualTo(2));
    }

    [Test]
    public void ChargeMismatchPenalty_AvoidsHighCost()
    {
        // Arrange: simulate charge mismatch with high penalty
        double[,] cost = new double[,]
        {
            { 0, 1000 },
            { 1000, 0 }
        };

        // Act
        int[] assignment = HungarianAlgorithm.FindAssignments(cost);

        // Assert
        Assert.That(assignment[0], Is.EqualTo(0));
        Assert.That(assignment[1], Is.EqualTo(1));
    }
}
