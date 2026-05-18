using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using EngineLayer;
using GuiFunctions;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using Assert = NUnit.Framework.Assert;

namespace Test.GuiTests;

[TestFixture]
public class MzLibExtensionsTests
{
    public PsmFromTsv TestPsm;
    public OsmFromTsv TestOsm;

    [OneTimeSetUp]
    public void OneTimeSetUp()
    {
        string psmPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "SequenceCoverageTestPSM.psmtsv");
        string osmPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "AllOligos.osmtsv");

        TestPsm = SpectrumMatchTsvReader.ReadTsv<PsmFromTsv>(psmPath, out _)[0];
        TestOsm = SpectrumMatchTsvReader.ReadTsv<OsmFromTsv>(osmPath, out _)[0];
    }


    [Test]
    public void ToViewModel_WithClassicDeconvolutionParameters_ReturnsClassicDeconParamsViewModel()
    {
        // Arrange
        var classicParams = new ClassicDeconvolutionParameters(1, 12, 10, 1);

        // Act
        var result = classicParams.ToViewModel();

        // Assert
        Assert.That(result, Is.InstanceOf<ClassicDeconParamsViewModel>());
    }

    [Test]
    public void ToViewModel_WithIsoDecDeconvolutionParameters_ReturnsClassicDeconParamsViewModel()
    {
        // Arrange
        var isoDecParams = new IsoDecDeconvolutionParameters();

        // Act
        var result = isoDecParams.ToViewModel();

        // Assert
        Assert.That(result, Is.InstanceOf<IsoDecDeconParamsViewModel>());
    }

    [Test]
    public void ToViewModel_WithUnsupportedParameters_ThrowsNotImplementedException()
    {
        // Arrange
        var unsupportedParams = new ExampleNewDeconvolutionParametersTemplate(1, 20);

        // Act & Assert
        Assert.That(() => unsupportedParams.ToViewModel(), Throws.TypeOf<NotImplementedException>());
    }

    [Test]
    public void MajorityWithin_ReturnsTrue_WhenMajorityWithinRange()
    {
        var range = new MzRange(10, 20);
        var values = new List<double> { 12, 15, 18, 25, 30 };
        // 3 out of 5 are within [10, 20]
        Assert.That(range.MajorityWithin(values), Is.True);
    }

    [Test]
    public void MajorityWithin_ReturnsFalse_WhenMajorityNotWithinRange()
    {
        var range = new MzRange(10, 20);
        var values = new List<double> { 5, 8, 12, 25, 30 };
        // Only 1 out of 5 is within [10, 20]
        Assert.That(range.MajorityWithin(values), Is.False);
    }
    [Test]
    public void IsPeptide_ReturnsTrue_ForPsmFromTsv()
    {
        // Arrange
        var psm = TestPsm;

        // Act
        var result = psm.IsPeptide();

        // Assert
        Assert.That(result, Is.True);
    }

    [Test]
    public void IsPeptide_ReturnsFalse_ForOsmFromTsv()
    {
        // Arrange
        var osm = TestOsm;

        // Act
        var result = osm.IsPeptide();

        // Assert
        Assert.That(result, Is.False);
    }

    [Test]
    public void GetDigestionProductLabel_ReturnsPeptide_ForPsmFromTsv()
    {
        // Arrange
        var psm = TestPsm;

        // Act
        var label = psm.GetDigestionProductLabel();

        // Assert
        Assert.That(label, Is.EqualTo("Peptide"));
    }

    [Test]
    public void GetDigestionProductLabel_ReturnsOligo_ForOsmFromTsv()
    {
        // Arrange
        var osm = TestOsm;

        // Act
        var label = osm.GetDigestionProductLabel();

        // Assert
        Assert.That(label, Is.EqualTo("Oligo"));
    }



    [Test]
    public void IsCrossLinkedPeptide_ReturnsTrue_ForCrossLinkedPsm()
    {
        var crossLinkedPsm = new PsmFromTsv(TestPsm, TestPsm.FullSequence, baseSequence: TestPsm.BaseSeq);

        // Attempt to set the property via reflection.
        var backingField = typeof(PsmFromTsv).GetField("<BetaPeptideBaseSequence>k__BackingField", BindingFlags.Instance | BindingFlags.NonPublic);
        if (backingField != null)
        {
            backingField.SetValue(crossLinkedPsm, "BETASEQ");
        }
        else
        {
            Assert.Fail("Failed to set BetaPeptideBaseSequence via reflection: neither a writable property nor the backing field was found.");
        }

        Assert.That(crossLinkedPsm.IsCrossLinkedPeptide(), Is.True);
    }

    [Test]
    public void IsCrossLinkedPeptide_ReturnsFalse_ForNonCrossLinkedPsm()
    {
        Assert.That(TestPsm.IsCrossLinkedPeptide(), Is.False);
    }

    [Test]
    public void ToBioPolymerWithSetMods_ReturnsPeptideWithSetMods_WhenIsPeptide()
    {
        var result = TestPsm.ToBioPolymerWithSetMods();
        Assert.That(result.GetType().Name, Is.EqualTo("PeptideWithSetModifications"));
    }

    [Test]
    public void ToBioPolymerWithSetMods_ReturnsOligoWithSetMods_WhenIsOligo()
    {
        var result = TestOsm.ToBioPolymerWithSetMods();
        Assert.That(result.GetType().Name, Is.EqualTo("OligoWithSetMods"));
    }

    [Test]
    public void ToBioPolymerWithSetMods_ParsesStartAndEndResidues()
    {
        // Use reflection to set StartAndEndResiduesInParentSequence
        var psm = TestPsm;
        var prop = typeof(SpectrumMatchFromTsv).GetProperty("StartAndEndResiduesInParentSequence");
        prop.SetValue(psm, "[2 to 5]");
        var result = psm.ToBioPolymerWithSetMods();
        // Use reflection to get OneBasedStartResidueInProtein/OneBasedEndResidueInProtein if available
        var startProp = result.GetType().GetProperty("OneBasedStartResidueInProtein");
        var endProp = result.GetType().GetProperty("OneBasedEndResidueInProtein");
        if (startProp != null && endProp != null)
        {
            Assert.That((int)startProp.GetValue(result), Is.EqualTo(2));
            Assert.That((int)endProp.GetValue(result), Is.EqualTo(5));
        }
    }

    [Test]
    public void ReplaceFullSequence_ReplacesForPeptide()
    {
        var replaced = TestPsm.ReplaceFullSequence("NEWSEQ", "NEWBASE");
        Assert.That(replaced, Is.InstanceOf<PsmFromTsv>());
        Assert.That(replaced.FullSequence, Is.EqualTo("NEWSEQ"));
    }

    [Test]
    public void ReplaceFullSequence_ReplacesForOligo()
    {
        var replaced = TestOsm.ReplaceFullSequence("NEWSEQ", "NEWBASE");
        Assert.That(replaced, Is.InstanceOf<OsmFromTsv>());
        Assert.That(replaced.FullSequence, Is.EqualTo("NEWSEQ"));
    }

    [Test]
    public void ProductsFromDissociationType_ReturnsPeptideDict_ForPeptide()
    {
        var dict = TestPsm.ProductsFromDissociationType();
        Assert.That(dict, Is.EqualTo(Omics.Fragmentation.Peptide.DissociationTypeCollection.ProductsFromDissociationType));
    }

    [Test]
    public void ProductsFromDissociationType_ReturnsOligoDict_ForOligo()
    {
        var dict = TestOsm.ProductsFromDissociationType();
        Assert.That(dict, Is.EqualTo(Omics.Fragmentation.Oligo.DissociationTypeCollection.ProductsFromDissociationType));
    }

    [Test]
    public void AllModsKnownModificationDictionary_ReturnsPeptideDict_ForPeptide()
    {
        var dict = TestPsm.AllModsKnownModificationDictionary();
        Assert.That(dict, Is.Not.Null);
        Assert.That(dict, Is.EqualTo(GlobalVariables.AllModsKnownDictionary));
    }

    [Test]
    public void AllModsKnownModificationDictionary_ReturnsOligoDict_ForOligo()
    {
        var dict = TestOsm.AllModsKnownModificationDictionary();
        Assert.That(dict, Is.Not.Null);
        Assert.That(dict, Is.EqualTo(GlobalVariables.AllRnaModsKnownDictionary));
    }

    [Test]
    public void GetStartAndEndPosition_ParsesMultipleRanges()
    {
        // Use reflection to set StartAndEndResiduesInParentSequence
        var psm = TestPsm;
        var prop = typeof(SpectrumMatchFromTsv).GetProperty("StartAndEndResiduesInParentSequence");
        prop.SetValue(psm, "[1to5]|[7to9]");

        var expected = new List<(int, int)> { (1, 5), (7, 9) };
        var positions = psm.GetStartAndEndPosition().ToList();
        Assert.That(positions.Count, Is.EqualTo(expected.Count));

        foreach (var pos in positions)
        {
            Assert.That(expected.Contains(pos));
        }
    }
}