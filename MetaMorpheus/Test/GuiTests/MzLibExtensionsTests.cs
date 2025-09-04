using System;
using System.Collections.Generic;
using System.IO;
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
        string osmPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData" , "AllOligos.osmtsv");

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
        Assert.That(label, Is.EqualTo("Pepide"));
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
}