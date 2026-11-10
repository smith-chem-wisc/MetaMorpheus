using System.Collections.Generic;
using GuiFunctions;
using GuiFunctions.MetaDraw;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.MetaDraw;

/// <summary>
/// The ChimeraGroupViewModelTests that build their group through the WPF ChimeraAnalysisTabViewModel.
/// Split out so the rest of the fixture compiles into GuiFunctions.Core.Tests and runs off Windows.
/// </summary>
public partial class ChimeraGroupViewModelTests
{
    [OneTimeSetUp]
    public static void SuppressMessageBoxes()
    {
        MessageBoxHelper.SuppressMessageBoxes = true;
    }

    [Test]
    public void ChimeraGroupViewModel_AssignIonColors_UsesLetterOption()
    {
        // Arrange
        var allPsms = ChimeraGroupViewModelTests.AllMatches;
        var dataFiles = new Dictionary<string, MsDataFile> { { "FXN3_tr1_032017-calib", ChimeraGroupViewModelTests.DataFile } };
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);
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
        var vm = new ChimeraAnalysisTabViewModel();
        vm.ProcessChimeraData(allPsms, dataFiles);
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

        var vm = new ChimeraAnalysisTabViewModel()
        {
            UseLetterOnly = true
        };
        vm.ProcessChimeraData(allPsms, dataFiles);

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
}
