using Chemistry;
using EngineLayer;
using GuiFunctions.Models;
using GuiFunctions.Util;
using NUnit.Framework;
using System;
using System.IO;
using System.Linq;

namespace Test.GuiTests;

[TestFixture]
public class CustomMIonLossTests
{
    private string _testFilePath;

    [SetUp]
    public void Setup()
    {
        // Use a temporary test file path
        _testFilePath = Path.Combine(Path.GetTempPath(), "TestCustomMIonLosses.txt");
        
        // Clear any existing cache and test file
        CustomMIonLossManager.ClearCache();
        if (File.Exists(_testFilePath))
        {
            File.Delete(_testFilePath);
        }
    }

    [TearDown]
    public void TearDown()
    {
        // Clean up test file
        if (File.Exists(_testFilePath))
        {
            File.Delete(_testFilePath);
        }
        CustomMIonLossManager.ClearCache();
    }

    #region CustomMIonLoss Tests

    [Test]
    public void CustomMIonLoss_Constructor_CreatesValidObject()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("H2O1");
        
        // Act
        var loss = new CustomMIonLoss("Water Loss", "-H2O", formula, AnalyteType.Peptide);
        
        // Assert
        Assert.That(loss.Name, Is.EqualTo("Water Loss"));
        Assert.That(loss.Annotation, Is.EqualTo("-H2O"));
        Assert.That(loss.ApplicableAnalyteType, Is.EqualTo(AnalyteType.Peptide));
        Assert.That(loss.MonoisotopicMass, Is.EqualTo(formula.MonoisotopicMass).Within(0.0001));
    }

    [Test]
    public void CustomMIonLoss_ToMIonLoss_ConvertsCorrectly()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("H3N1");
        var customLoss = new CustomMIonLoss("Ammonia Loss", "-NH3", formula, AnalyteType.Peptide);
        
        // Act
        var mIonLoss = customLoss.ToMIonLoss();
        
        // Assert
        Assert.That(mIonLoss.Name, Is.EqualTo("Ammonia Loss"));
        Assert.That(mIonLoss.Annotation, Is.EqualTo("-NH3"));
        Assert.That(mIonLoss.MonoisotopicMass, Is.EqualTo(customLoss.MonoisotopicMass).Within(0.0001));
    }

    [Test]
    public void CustomMIonLoss_FromMIonLoss_CreatesCorrectly()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("H1P1O3");
        var mIonLoss = new Omics.Fragmentation.MIonLoss("Phosphate Loss", "-P", formula);
        
        // Act
        var customLoss = CustomMIonLoss.FromMIonLoss(mIonLoss, AnalyteType.Oligo);
        
        // Assert
        Assert.That(customLoss.Name, Is.EqualTo(mIonLoss.Name));
        Assert.That(customLoss.Annotation, Is.EqualTo(mIonLoss.Annotation));
        Assert.That(customLoss.ApplicableAnalyteType, Is.EqualTo(AnalyteType.Oligo));
        Assert.That(customLoss.MonoisotopicMass, Is.EqualTo(mIonLoss.MonoisotopicMass).Within(0.0001));
    }

    [Test]
    public void CustomMIonLoss_Serialize_CreatesCorrectFormat()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("H2O1");
        var loss = new CustomMIonLoss("Water Loss", "-H2O", formula, AnalyteType.Peptide);
        
        // Act
        var serialized = loss.Serialize();
        
        // Assert
        Assert.That(serialized, Is.EqualTo("Water Loss|-H2O|H2O|Peptide"));
    }

    [Test]
    public void CustomMIonLoss_Deserialize_ParsesCorrectly()
    {
        // Arrange
        var line = "Water Loss|-H2O|H2O1|Peptide";
        
        // Act
        var loss = CustomMIonLoss.Deserialize(line);
        
        // Assert
        Assert.That(loss.Name, Is.EqualTo("Water Loss"));
        Assert.That(loss.Annotation, Is.EqualTo("-H2O"));
        Assert.That(loss.ApplicableAnalyteType, Is.EqualTo(AnalyteType.Peptide));
        Assert.That(loss.MonoisotopicMass, Is.GreaterThan(0));
    }

    [Test]
    public void CustomMIonLoss_Deserialize_InvalidFormat_ThrowsException()
    {
        // Arrange
        var invalidLine = "Water Loss|-H2O|H2O1"; // Missing AnalyteType
        
        // Act & Assert
        Assert.Throws<FormatException>(() => CustomMIonLoss.Deserialize(invalidLine));
    }

    [Test]
    public void CustomMIonLoss_Deserialize_InvalidAnalyteType_ThrowsException()
    {
        // Arrange
        var invalidLine = "Water Loss|-H2O|H2O1|InvalidType";
        
        // Act & Assert
        Assert.Throws<FormatException>(() => CustomMIonLoss.Deserialize(invalidLine));
    }

    [Test]
    public void CustomMIonLoss_Deserialize_InvalidFormula_ThrowsException()
    {
        // Arrange
        var invalidLine = "Water Loss|-H2O|InvalidFormula123|Peptide";
        
        // Act & Assert
        Assert.Throws<FormatException>(() => CustomMIonLoss.Deserialize(invalidLine));
    }

    [Test]
    public void CustomMIonLoss_IsApplicableToCurrentMode_PeptideMode_ReturnsTrue()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("H2O1");
        var loss = new CustomMIonLoss("Water Loss", "-H2O", formula, AnalyteType.Peptide);
        
        // Act
        var isApplicable = loss.IsApplicableToCurrentMode(AnalyteType.Peptide);
        
        // Assert
        Assert.That(isApplicable, Is.True);
    }

    [Test]
    public void CustomMIonLoss_IsApplicableToCurrentMode_OligoMode_ReturnsTrue()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("H1P1O3");
        var loss = new CustomMIonLoss("Phosphate Loss", "-P", formula, AnalyteType.Oligo);
        
        // Act
        var isApplicable = loss.IsApplicableToCurrentMode(AnalyteType.Oligo);
        
        // Assert
        Assert.That(isApplicable, Is.True);
    }

    [Test]
    public void CustomMIonLoss_IsApplicableToCurrentMode_WrongMode_ReturnsFalse()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("H2O1");
        var loss = new CustomMIonLoss("Water Loss", "-H2O", formula, AnalyteType.Peptide);
        
        // Act
        var isApplicable = loss.IsApplicableToCurrentMode(AnalyteType.Oligo);
        
        // Assert
        Assert.That(isApplicable, Is.False);
    }

    [Test]
    public void CustomMIonLoss_SerializeDeserialize_RoundTrip_Successful()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("C2H4O2");
        var originalLoss = new CustomMIonLoss("Acetic Acid Loss", "-CH3COOH", formula, AnalyteType.Peptide);
        
        // Act
        var serialized = originalLoss.Serialize();
        var deserializedLoss = CustomMIonLoss.Deserialize(serialized);
        
        // Assert
        Assert.That(deserializedLoss.Name, Is.EqualTo(originalLoss.Name));
        Assert.That(deserializedLoss.Annotation, Is.EqualTo(originalLoss.Annotation));
        Assert.That(deserializedLoss.ApplicableAnalyteType, Is.EqualTo(originalLoss.ApplicableAnalyteType));
        Assert.That(deserializedLoss.MonoisotopicMass, Is.EqualTo(originalLoss.MonoisotopicMass).Within(0.0001));
    }

    #endregion

    #region CustomMIonLossManager Tests

    [Test]
    public void CustomMIonLossManager_GetCustomMIonLossFilePath_ReturnsValidPath()
    {
        // Act
        var filePath = CustomMIonLossManager.GetCustomMIonLossFilePath();
        
        // Assert
        Assert.That(filePath, Is.Not.Null);
        Assert.That(filePath, Does.Contain("Mods"));
        Assert.That(filePath, Does.EndWith("CustomMIonLosses.txt"));
    }

    [Test]
    public void CustomMIonLossManager_AddCustomMIonLoss_AddsSuccessfully()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("CO2");
        var newLoss = new CustomMIonLoss("CO2 Loss", "-CO2", formula, AnalyteType.Peptide);
        
        // Act
        CustomMIonLossManager.AddCustomMIonLoss(newLoss);
        var losses = CustomMIonLossManager.LoadCustomMIonLosses();
        
        // Assert
        Assert.That(losses.Any(l => l.Name == "CO2 Loss"), Is.True);
        Assert.That(losses.Any(l => l.Annotation == "-CO2"), Is.True);
    }

    [Test]
    public void CustomMIonLossManager_AddCustomMIonLoss_Duplicate_ThrowsException()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("H2O1");
        var loss1 = new CustomMIonLoss("Test Loss", "-TEST", formula, AnalyteType.Peptide);
        var loss2 = new CustomMIonLoss("Test Loss", "-TEST2", formula, AnalyteType.Peptide);
        
        // Act
        CustomMIonLossManager.AddCustomMIonLoss(loss1);
        
        // Assert
        Assert.Throws<InvalidOperationException>(() => CustomMIonLossManager.AddCustomMIonLoss(loss2));
    }

    [Test]
    public void CustomMIonLossManager_AddCustomMIonLoss_SameNameDifferentMode_Succeeds()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("H2O1");
        var peptideLoss = new CustomMIonLoss("Multi Mode Loss", "-MMM", formula, AnalyteType.Peptide);
        var oligoLoss = new CustomMIonLoss("Multi Mode Loss", "-MMM", formula, AnalyteType.Oligo);
        
        // Act
        CustomMIonLossManager.AddCustomMIonLoss(peptideLoss);
        CustomMIonLossManager.AddCustomMIonLoss(oligoLoss);
        var losses = CustomMIonLossManager.LoadCustomMIonLosses();
        
        // Assert
        var multiModeLosses = losses.Where(l => l.Name == "Multi Mode Loss").ToList();
        Assert.That(multiModeLosses.Count, Is.EqualTo(2));
        Assert.That(multiModeLosses.Any(l => l.ApplicableAnalyteType == AnalyteType.Peptide), Is.True);
        Assert.That(multiModeLosses.Any(l => l.ApplicableAnalyteType == AnalyteType.Oligo), Is.True);
    }

    [Test]
    public void CustomMIonLossManager_SaveAndLoad_PreservesData()
    {
        // Arrange
        var formula1 = ChemicalFormula.ParseFormula("H2O1");
        var formula2 = ChemicalFormula.ParseFormula("H3N1");
        var losses = new[]
        {
            new CustomMIonLoss("Water Loss", "-H2O", formula1, AnalyteType.Peptide),
            new CustomMIonLoss("Ammonia Loss", "-NH3", formula2, AnalyteType.Oligo)
        };
        
        // Act
        CustomMIonLossManager.SaveCustomMIonLosses(losses);
        CustomMIonLossManager.ClearCache(); // Force reload from file
        var loadedLosses = CustomMIonLossManager.LoadCustomMIonLosses();
        
        // Assert
        Assert.That(loadedLosses.Count, Is.GreaterThanOrEqualTo(2));
        Assert.That(loadedLosses.Any(l => l.Name == "Water Loss" && l.ApplicableAnalyteType == AnalyteType.Peptide), Is.True);
        Assert.That(loadedLosses.Any(l => l.Name == "Ammonia Loss" && l.ApplicableAnalyteType == AnalyteType.Oligo), Is.True);
    }

    [Test]
    public void CustomMIonLossManager_ClearCache_ForcesReload()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("CO1");
        var newLoss = new CustomMIonLoss("CO Loss", "-CO", formula, AnalyteType.Peptide);
        
        // Act
        var losses1 = CustomMIonLossManager.LoadCustomMIonLosses();
        var initialCount = losses1.Count;
        
        CustomMIonLossManager.AddCustomMIonLoss(newLoss);
        CustomMIonLossManager.ClearCache();
        
        var losses2 = CustomMIonLossManager.LoadCustomMIonLosses();
        
        // Assert
        Assert.That(losses2.Count, Is.EqualTo(initialCount + 1));
        Assert.That(losses2.Any(l => l.Name == "CO Loss"), Is.True);
    }

    [Test]
    public void CustomMIonLossManager_DeleteCustomMIonLoss_RemovesSuccessfully()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("SO3");
        var lossToDelete = new CustomMIonLoss("SO3 Loss", "-SO3", formula, AnalyteType.Peptide);
        CustomMIonLossManager.AddCustomMIonLoss(lossToDelete);
        
        // Act
        CustomMIonLossManager.DeleteCustomMIonLoss("SO3 Loss", AnalyteType.Peptide);
        var losses = CustomMIonLossManager.LoadCustomMIonLosses();
        
        // Assert
        Assert.That(losses.Any(l => l.Name == "SO3 Loss" && l.ApplicableAnalyteType == AnalyteType.Peptide), Is.False);
    }

    [Test]
    public void CustomMIonLossManager_DeleteCustomMIonLoss_NonExistent_DoesNotThrow()
    {
        // Act & Assert
        Assert.DoesNotThrow(() => CustomMIonLossManager.DeleteCustomMIonLoss("NonExistent Loss", AnalyteType.Peptide));
    }

    [Test]
    public void CustomMIonLossManager_GetAllMIonLossesForAnalyteType_Peptide_ReturnsCorrectLosses()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("H2O1");
        var peptideLoss = new CustomMIonLoss("Peptide Only Loss", "-PEP", formula, AnalyteType.Peptide);
        CustomMIonLossManager.AddCustomMIonLoss(peptideLoss);
        
        // Act
        var peptideLosses = CustomMIonLossManager.GetAllMIonLossesForAnalyteType(AnalyteType.Peptide);
        
        // Assert
        Assert.That(peptideLosses.Any(l => l.Name == "Peptide Only Loss"), Is.True);
    }

    [Test]
    public void CustomMIonLossManager_GetAllMIonLossesForAnalyteType_Oligo_ReturnsCorrectLosses()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("H1P1O3");
        var oligoLoss = new CustomMIonLoss("Oligo Only Loss", "-OLI", formula, AnalyteType.Oligo);
        CustomMIonLossManager.AddCustomMIonLoss(oligoLoss);
        
        // Act
        var oligoLosses = CustomMIonLossManager.GetAllMIonLossesForAnalyteType(AnalyteType.Oligo);
        
        // Assert
        Assert.That(oligoLosses.Any(l => l.Name == "Oligo Only Loss"), Is.True);
    }

    [Test]
    public void CustomMIonLossManager_GetAllMIonLossesForAnalyteType_FiltersCorrectly()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("H2O1");
        var peptideLoss = new CustomMIonLoss("Peptide Specific", "-PEP", formula, AnalyteType.Peptide);
        var oligoLoss = new CustomMIonLoss("Oligo Specific", "-OLI", formula, AnalyteType.Oligo);
        
        CustomMIonLossManager.AddCustomMIonLoss(peptideLoss);
        CustomMIonLossManager.AddCustomMIonLoss(oligoLoss);
        
        // Act
        var peptideLosses = CustomMIonLossManager.GetAllMIonLossesForAnalyteType(AnalyteType.Peptide);
        var oligoLosses = CustomMIonLossManager.GetAllMIonLossesForAnalyteType(AnalyteType.Oligo);
        
        // Assert
        Assert.That(peptideLosses.Any(l => l.Name == "Peptide Specific"), Is.True);
        Assert.That(peptideLosses.Any(l => l.Name == "Oligo Specific"), Is.False);
        
        Assert.That(oligoLosses.Any(l => l.Name == "Oligo Specific"), Is.True);
        Assert.That(oligoLosses.Any(l => l.Name == "Peptide Specific"), Is.False);
    }

    [Test]
    public void CustomMIonLossManager_LoadCustomMIonLosses_SkipsComments()
    {
        // Arrange
        var filePath = CustomMIonLossManager.GetCustomMIonLossFilePath();
        var lines = new[]
        {
            "# This is a comment",
            "Water Loss|-H2O|H2O1|Peptide",
            "# Another comment",
            "",
            "Ammonia Loss|-NH3|H3N1|Peptide"
        };
        File.WriteAllLines(filePath, lines);
        CustomMIonLossManager.ClearCache();
        
        // Act
        var losses = CustomMIonLossManager.LoadCustomMIonLosses();
        
        // Assert
        Assert.That(losses.Count, Is.EqualTo(2));
        Assert.That(losses.All(l => !l.Name.StartsWith("#")), Is.True);
    }

    [Test]
    public void CustomMIonLossManager_LoadCustomMIonLosses_UsesCaching()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("H2O1");
        var loss = new CustomMIonLoss("Cache Test", "-CACHE", formula, AnalyteType.Peptide);
        CustomMIonLossManager.AddCustomMIonLoss(loss);
        
        // Act
        var losses1 = CustomMIonLossManager.LoadCustomMIonLosses();
        var losses2 = CustomMIonLossManager.LoadCustomMIonLosses(); // Should use cache
        
        // Assert
        Assert.That(losses1.Count, Is.EqualTo(losses2.Count));
        Assert.That(losses1.Select(l => l.Name), Is.EquivalentTo(losses2.Select(l => l.Name)));
    }

    [Test]
    public void CustomMIonLossManager_LoadCustomMIonLosses_InvalidLine_SkipsAndContinues()
    {
        // Arrange
        var filePath = CustomMIonLossManager.GetCustomMIonLossFilePath();
        var lines = new[]
        {
            "Water Loss|-H2O|H2O1|Peptide",
            "Invalid|Line", // Invalid format
            "Ammonia Loss|-NH3|H3N1|Peptide"
        };
        File.WriteAllLines(filePath, lines);
        CustomMIonLossManager.ClearCache();
        
        // Act
        var losses = CustomMIonLossManager.LoadCustomMIonLosses();
        
        // Assert
        Assert.That(losses.Count, Is.EqualTo(2)); // Should skip invalid line
        Assert.That(losses.Any(l => l.Name == "Water Loss"), Is.True);
        Assert.That(losses.Any(l => l.Name == "Ammonia Loss"), Is.True);
    }

    #endregion
}
