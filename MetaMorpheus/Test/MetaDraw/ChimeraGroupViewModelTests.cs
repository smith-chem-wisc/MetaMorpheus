using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Chemistry;
using EngineLayer;
using GuiFunctions;
using GuiFunctions.MetaDraw;
using MassSpectrometry;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using Omics.Fragmentation;
using OxyPlot;
using Readers;
using Path = System.IO.Path;

namespace Test.MetaDraw;

[TestFixture]
[ExcludeFromCodeCoverage]
public class ChimeraGroupViewModelTests
{
    public static ChimeraTestCase OneProteinTwoProteoformChimeraGroup;
    public static ChimeraTestCase TwoProteinsTwoProteoformChimeraGroup;
    public static List<SpectrumMatchFromTsv> AllMatches;
    public static List<SpectrumMatchFromTsv> AllMatchesMutable;
    public static MsDataFile DataFile;
    public static string TestExportDirectory => Path.Combine(TestContext.CurrentContext.TestDirectory, "MetaDraw", "ChimeraPlottingTests");
    public record ChimeraTestCase(ChimeraGroupViewModel ChimeraGroup, Dictionary<OxyColor, List<MatchedFragmentIon>> ExpectedIonsByColor);

    static ChimeraGroupViewModelTests()
    {
        var psmPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TDGPTMDSearchResults.psmtsv");
        string msDataPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TDGPTMDSearchSingleSpectra.mzML");

        DataFile = MsDataFileReader.GetDataFile(msDataPath).LoadAllStaticData();
        DataFile.InitiateDynamicConnection();
        AllMatches = SpectrumMatchTsvReader.ReadTsv(psmPath, out var warnings).Where(p => p.Ms2ScanNumber <= DataFile.Scans.Length).ToList();
        AllMatchesMutable = SpectrumMatchTsvReader.ReadTsv(psmPath, out warnings).Where(p => p.Ms2ScanNumber <= DataFile.Scans.Length).ToList();
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
        MessageBoxHelper.SuppressMessageBoxes = true;
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
        // Clean up the export directory after tests
        if (Directory.Exists(TestExportDirectory))
        {
            Directory.Delete(TestExportDirectory, true);
        }
    }

    [Test]
    public static void FragmentCounts_OneProteinTwoProteoforms()
    {
        // Arrange
        var chimeraGroup = OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        Assert.That(chimeraGroup.Count, Is.EqualTo(2));
        Assert.That(chimeraGroup.ProteinCount, Is.EqualTo(1));
        Assert.That(chimeraGroup.OneBasedPrecursorScanNumber, Is.EqualTo(chimeraGroup.ChimericPsms.First().Psm.PrecursorScanNum));
        Assert.That(chimeraGroup.PrimarySequenceCount, Is.EqualTo(1));

        var terminalFragments = chimeraGroup.ChimericPsms.SelectMany(p => p.Psm.MatchedIons)
            .Where(p => !p.IsInternalFragment)
            .ToList();
        Assert.That(terminalFragments.Count, Is.EqualTo(chimeraGroup.TotalFragments), chimeraGroup.TotalFragments + " terminal fragments should match total fragments in the group.");

        var uniqueFragments = terminalFragments.Select(p => (p.NeutralTheoreticalProduct.Annotation, p.Mz.RoundedDouble(1)))
            .Distinct()
            .ToList();
        Assert.That(uniqueFragments.Count, Is.EqualTo(chimeraGroup.UniqueFragments), "Unique fragments should match total unique fragments in the group.");
    }

    [Test]
    public static void FragmentCounts_TwoProteinsTwoProteoforms()
    {
        // Arrange
        var chimeraGroup = TwoProteinsTwoProteoformChimeraGroup.ChimeraGroup;
        Assert.That(chimeraGroup.Count, Is.EqualTo(2));
        Assert.That(chimeraGroup.ProteinCount, Is.EqualTo(2));
        Assert.That(chimeraGroup.PrimarySequenceCount, Is.EqualTo(2));

        var terminalFragments = chimeraGroup.ChimericPsms.SelectMany(p => p.Psm.MatchedIons)
            .Where(p => !p.IsInternalFragment)
            .ToList();
        Assert.That(terminalFragments.Count, Is.EqualTo(chimeraGroup.TotalFragments), chimeraGroup.TotalFragments + " terminal fragments should match total fragments in the group.");

        var uniqueFragments = terminalFragments.Select(p => (p.NeutralTheoreticalProduct.Annotation, p.Mz.RoundedDouble(1)))
            .Distinct()
            .ToList();
        Assert.That(uniqueFragments.Count, Is.EqualTo(chimeraGroup.UniqueFragments), "Unique fragments should match total unique fragments in the group.");
    }

    [Test]
    public static void PrecursorAssignmentIsCorrect_OneProteinTwoProteoforms()
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
    public static void PrecursorAssignmentIsCorrect_TwoProteinsTwoProteoforms()
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
    public static void PrecursorColorAssignmentIsCorrect_OneProteinTwoProteoforms()
    {
        // Arrange
        var chimeraGroup = OneProteinTwoProteoformChimeraGroup.ChimeraGroup;
        Assert.That(chimeraGroup.Count, Is.EqualTo(2));

        var envelope1 = chimeraGroup.ChimericPsms[0].PrecursorEnvelope;
        var envelope2 = chimeraGroup.ChimericPsms[1].PrecursorEnvelope;

        var envelope1Peaks = envelope1.Peaks;
        var envelope2Peaks = envelope2.Peaks;

        var sharedPeakColor = ChimeraGroupViewModel.ColorByProteinDictionary[0][0];
        var firstColor = ChimeraGroupViewModel.ColorByProteinDictionary[0][1];
        var secondColor = ChimeraGroupViewModel.ColorByProteinDictionary[0][2];

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
    public static void PrecursorColorAssignmentIsCorrect_TwoProteinsTwoProteoforms()
    {
        // Arrange
        var chimeraGroup = TwoProteinsTwoProteoformChimeraGroup.ChimeraGroup;
        Assert.That(chimeraGroup.Count, Is.EqualTo(2));

        var envelope1 = chimeraGroup.ChimericPsms[0].PrecursorEnvelope;
        var envelope2 = chimeraGroup.ChimericPsms[1].PrecursorEnvelope;

        var envelope1Peaks = envelope1.Peaks;
        var envelope2Peaks = envelope2.Peaks;

        var sharedPeakColor = OxyColors.Black;
        var firstProteinColor = ChimeraGroupViewModel.ColorByProteinDictionary[0][0];
        var firstProteoformOfFirstProteinColor = ChimeraGroupViewModel.ColorByProteinDictionary[0][1];
        var secondProteinColor = ChimeraGroupViewModel.ColorByProteinDictionary[1][0];
        var firstProteoformOfSecondProteinColor = ChimeraGroupViewModel.ColorByProteinDictionary[1][1];

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
    public static void FragmentColorAssignmentIsCorrect_OneProteinTwoProteoforms()
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
    public static void FragmentColorAssignmentIsCorrect_TwoProteinsTwoProteoforms()
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

    [Test]
    public void ChimeraGroupViewModel_AssignIonColors_UsesLetterOption()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles);
        var group = vm.ChimeraGroupViewModels[0];

        // Act
        group.AssignIonColors(useLetterOnly: true);

        // Assert
        foreach (var kvp in group.PrecursorIonsByColor)
        {
            foreach (var tuple in kvp.Value)
            {
                // If annotation is not "Miso", it should contain a letter from the Letters list
                if (!tuple.Item2.Contains("Miso"))
                {
                    Assert.That(tuple.Item2, Is.Not.Empty, "Annotation should not be empty when using letter only.");
                    Assert.That(tuple.Item2.Length, Is.EqualTo(1), "Annotation should be a single letter when using letter only.");
                    Assert.That("ABCDEFGHIJKLMNOPQRSTUVWXYZ".Contains(tuple.Item2), Is.True, "Annotation should be a letter.");
                }
            }
        }
    }

    [Test]
    public void AssignIonColors_ProduceExpectedAnnotations()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles);
        var group = vm.ChimeraGroupViewModels[0];

        // Act
        group.AssignIonColors(useLetterOnly: false);

        // Assert
        foreach (var kvp in group.PrecursorIonsByColor)
        {
            foreach (var tuple in kvp.Value)
            {
                // If annotation is not "Miso", it should contain the expected details
                if (!tuple.Item2.Contains("Miso"))
                {
                    Assert.That(tuple.Item2, Does.Contain("Charge ="), "Annotation should contain charge info.");
                    Assert.That(tuple.Item2, Does.Contain("m/z ="), "Annotation should contain m/z info.");
                    Assert.That(tuple.Item2, Does.Contain("Mono Mass ="), "Annotation should contain mono mass info.");
                    Assert.That(tuple.Item2, Does.Contain("Protein ="), "Annotation should contain protein info.");
                }
            }
        }
    }

    [Test]
    public void AssignIonColors_ProduceLetterOnlyAnnotations()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel(allPsms, dataFiles)
        {
            UseLetterOnly = true
        };
        vm.UseLetterOnly = true;
        var group = vm.ChimeraGroupViewModels[0];

        // Act
        group.AssignIonColors(useLetterOnly: true);

        // Assert
        foreach (var kvp in group.PrecursorIonsByColor)
        {
            foreach (var tuple in kvp.Value)
            {
                // If annotation is not "Miso", it should be a single uppercase letter
                if (!tuple.Item2.Contains("Miso"))
                {
                    Assert.That(tuple.Item2, Is.Not.Empty, "Annotation should not be empty when using letter only.");
                    Assert.That(tuple.Item2.Length, Is.EqualTo(1), "Annotation should be a single letter when using letter only.");
                    Assert.That("ABCDEFGHIJKLMNOPQRSTUVWXYZ".Contains(tuple.Item2), Is.True, "Annotation should be a letter.");
                }
            }
        }
    }

    [Test]
    public void ChimeraGroupViewModel_IEnumerable_EnumeratesAllChimericPsms()
    {
        // Arrange
        var chimeraGroup = ChimeraGroupViewModelTests.OneProteinTwoProteoformChimeraGroup.ChimeraGroup;

        // Act
        var enumerated = chimeraGroup.ToList();

        // Assert
        Assert.That(enumerated.Count, Is.EqualTo(chimeraGroup.ChimericPsms.Count));
        CollectionAssert.AreEqual(chimeraGroup.ChimericPsms, enumerated);
    }

    [Test]
    public void ChimeraGroupViewModel_IEnumerable_NonGenericEnumeratesAllChimericPsms()
    {
        // Arrange
        var chimeraGroup = ChimeraGroupViewModelTests.TwoProteinsTwoProteoformChimeraGroup.ChimeraGroup;

        // Act
        var enumerated = new List<ChimericSpectralMatchModel>();
        foreach (object item in (System.Collections.IEnumerable)chimeraGroup)
        {
            enumerated.Add((ChimericSpectralMatchModel)item);
        }

        // Assert
        Assert.That(enumerated.Count, Is.EqualTo(chimeraGroup.ChimericPsms.Count));
        CollectionAssert.AreEqual(chimeraGroup.ChimericPsms, enumerated);
    }
}