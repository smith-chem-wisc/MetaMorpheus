using System;
using System.Collections.Generic;
using GuiFunctions;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Assert = NUnit.Framework.Assert;

namespace Test.GuiTests;

[TestFixture]
public class MzLibExtensionsTests
{
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
}