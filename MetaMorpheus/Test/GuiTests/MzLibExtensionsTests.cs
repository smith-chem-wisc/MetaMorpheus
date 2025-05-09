using System;
using GuiFunctions;
using MassSpectrometry;
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
}