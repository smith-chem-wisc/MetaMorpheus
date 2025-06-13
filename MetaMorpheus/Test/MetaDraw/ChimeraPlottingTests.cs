using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EngineLayer;
using GuiFunctions;
using MassSpectrometry;
using NUnit.Framework;
using OxyPlot.Wpf;
using Readers;

namespace Test.MetaDraw;

[TestFixture]
[ExcludeFromCodeCoverage]
public class ChimeraPlottingTests
{
    public static ChimeraGroupViewModel TestChimeraGroup;
    public static List<SpectrumMatchFromTsv> AllMatches;
    public static MsDataFile DataFile;
    public static string TestExportDirectory => Path.Combine(TestContext.CurrentContext.TestDirectory, "MetaDraw", "ChimeraPlottingTests");

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

        var psmPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TDGPTMDSearchResults.psmtsv");
        string msDataPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TDGPTMDSearchSingleSpectra.mzML");


        DataFile = MsDataFileReader.GetDataFile(msDataPath).LoadAllStaticData();
        DataFile.InitiateDynamicConnection();
        AllMatches = SpectrumMatchTsvReader.ReadTsv(psmPath, out var warnings).Where(p => p.Ms2ScanNumber <= DataFile.Scans.Length).ToList();
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
