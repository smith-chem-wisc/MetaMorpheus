using GuiFunctions;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.GuiTests;

[TestFixture]
public class ClassicDeconParamsViewModelTest
{
    private ClassicDeconParamsViewModel _viewModel;
    private ClassicDeconvolutionParameters _parameters;

    [SetUp]
    public void SetUp()
    {
        _parameters = new ClassicDeconvolutionParameters(1, 12, 10, 0.5);
        _viewModel = new ClassicDeconParamsViewModel(_parameters);
    }

    [Test]
    public void TestDeconvolutionTolerancePpm()
    {
        Assert.That(_viewModel.DeconvolutionTolerancePpm, Is.EqualTo(10.0));
        _viewModel.DeconvolutionTolerancePpm = 20.0;
        Assert.That(_viewModel.DeconvolutionTolerancePpm, Is.EqualTo(20.0));
        Assert.That(_parameters.DeconvolutionTolerancePpm, Is.EqualTo(20.0));
    }

    [Test]
    public void TestIntensityRatioLimit()
    {
        Assert.That(_viewModel.IntensityRatioLimit, Is.EqualTo(0.5));
        _viewModel.IntensityRatioLimit = 1.0;
        Assert.That(_viewModel.IntensityRatioLimit, Is.EqualTo(1.0));
        Assert.That(_parameters.IntensityRatioLimit, Is.EqualTo(1.0));
    }

    [Test]
    public void TestToString()
    {
        Assert.That(_viewModel.ToString(), Is.EqualTo("Classic"));
    }
}