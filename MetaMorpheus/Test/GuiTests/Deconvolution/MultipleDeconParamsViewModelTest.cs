using System.Linq;
using GuiFunctions;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.GuiTests.Deconvolution;

[TestFixture]
public class MultipleDeconParamsViewModelTest
{
    private MultipleDeconParamsViewModel _viewModel;
    private ClassicDeconvolutionParameters _innerParams;
    private MultipleDeconParameters _parameters;

    [SetUp]
    public void SetUp()
    {
        _innerParams = new ClassicDeconvolutionParameters(1, 12, 4, 3, Polarity.Positive);
        _parameters = new MultipleDeconParameters(
            [_innerParams],
            _innerParams.MinAssumedChargeState,
            _innerParams.MaxAssumedChargeState,
            _innerParams.Polarity,
            _innerParams.AverageResidueModel,
            _innerParams.ExpectedIsotopeSpacing);
        _viewModel = new MultipleDeconParamsViewModel(_parameters);
    }

    [Test]
    public void Constructor_CreatesWithOneSubParameter()
    {
        Assert.That(_viewModel.SubParameters.Count, Is.EqualTo(1));
    }

    [Test]
    public void Constructor_SubParameterIsClassic()
    {
        var sub = _viewModel.SubParameters.Single();
        Assert.That(sub, Is.InstanceOf<ClassicDeconParamsViewModel>());
        Assert.That(sub.DeconvolutionType, Is.EqualTo(DeconvolutionType.ClassicDeconvolution));
    }

    [Test]
    public void ToString_ReturnsMultiple()
    {
        Assert.That(_viewModel.ToString(), Is.EqualTo("Multiple"));
    }

    [Test]
    public void AvailableSubTypes_ContainsClassicAndIsoDec()
    {
        Assert.That(_viewModel.AvailableSubTypes, Does.Contain(DeconvolutionType.ClassicDeconvolution));
        Assert.That(_viewModel.AvailableSubTypes, Does.Contain(DeconvolutionType.IsoDecDeconvolution));
        Assert.That(_viewModel.AvailableSubTypes.Count, Is.EqualTo(2));
    }

    [Test]
    public void AddSubType_AddsClassicSubParameter()
    {
        _viewModel.AddSubType(DeconvolutionType.ClassicDeconvolution);

        Assert.That(_viewModel.SubParameters.Count, Is.EqualTo(2));
        Assert.That(_viewModel.SubParameters.Last().DeconvolutionType, Is.EqualTo(DeconvolutionType.ClassicDeconvolution));
    }

    [Test]
    public void AddSubType_AddsIsoDecSubParameter()
    {
        _viewModel.AddSubType(DeconvolutionType.IsoDecDeconvolution);

        Assert.That(_viewModel.SubParameters.Count, Is.EqualTo(2));
        Assert.That(_viewModel.SubParameters.Last().DeconvolutionType, Is.EqualTo(DeconvolutionType.IsoDecDeconvolution));
    }

    [Test]
    public void AddSubType_InheritsCurrentSharedValues()
    {
        _viewModel.SharedMinAssumedChargeState = 2;
        _viewModel.SharedMaxAssumedChargeState = 10;
        _viewModel.AddSubType(DeconvolutionType.IsoDecDeconvolution);

        var added = _viewModel.SubParameters.Last();
        Assert.That(added.MinAssumedChargeState, Is.EqualTo(2));
        Assert.That(added.MaxAssumedChargeState, Is.EqualTo(10));
        Assert.That(added.Polarity, Is.EqualTo(Polarity.Positive));
    }

    [Test]
    public void RemoveSubType_RemovesSubParameter()
    {
        _viewModel.AddSubType(DeconvolutionType.IsoDecDeconvolution);
        Assert.That(_viewModel.SubParameters.Count, Is.EqualTo(2));

        _viewModel.RemoveSubType(_viewModel.SubParameters.Last());

        Assert.That(_viewModel.SubParameters.Count, Is.EqualTo(1));
        Assert.That(_viewModel.SubParameters.Single().DeconvolutionType, Is.EqualTo(DeconvolutionType.ClassicDeconvolution));
    }

    [Test]
    public void RemoveSubType_BlockWhenOnlyOne()
    {
        Assert.That(_viewModel.SubParameters.Count, Is.EqualTo(1));

        _viewModel.RemoveSubType(_viewModel.SubParameters.Single());

        Assert.That(_viewModel.SubParameters.Count, Is.EqualTo(1));
    }

    [Test]
    public void RemoveSubType_RebuildsParameters()
    {
        _viewModel.AddSubType(DeconvolutionType.IsoDecDeconvolution);
        var originalParameters = _viewModel.Parameters;
        _viewModel.RemoveSubType(_viewModel.SubParameters.Last());

        Assert.That(_viewModel.Parameters, Is.Not.SameAs(originalParameters));
    }

    [Test]
    public void SharedMinAssumedChargeState_PropagatesToSubParameters()
    {
        _viewModel.AddSubType(DeconvolutionType.IsoDecDeconvolution);
        _viewModel.SharedMinAssumedChargeState = 2;

        foreach (var sub in _viewModel.SubParameters)
            Assert.That(sub.MinAssumedChargeState, Is.EqualTo(2));
        Assert.That(_viewModel.MinAssumedChargeState, Is.EqualTo(2));
    }

    [Test]
    public void SharedMaxAssumedChargeState_PropagatesToSubParameters()
    {
        _viewModel.AddSubType(DeconvolutionType.IsoDecDeconvolution);
        _viewModel.SharedMaxAssumedChargeState = 15;

        foreach (var sub in _viewModel.SubParameters)
            Assert.That(sub.MaxAssumedChargeState, Is.EqualTo(15));
        Assert.That(_viewModel.MaxAssumedChargeState, Is.EqualTo(15));
    }

    [Test]
    public void Polarity_PropagatesToSubParameters()
    {
        _viewModel.AddSubType(DeconvolutionType.IsoDecDeconvolution);
        _viewModel.Polarity = Polarity.Negative;

        foreach (var sub in _viewModel.SubParameters)
            Assert.That(sub.Polarity, Is.EqualTo(Polarity.Negative));
        Assert.That(_viewModel.Polarity, Is.EqualTo(Polarity.Negative));
    }

    [Test]
    public void SelectedAddType_DefaultsToClassic()
    {
        Assert.That(_viewModel.SelectedAddType, Is.EqualTo(DeconvolutionType.ClassicDeconvolution));
    }

    [Test]
    public void Parameters_ReturnsMultipleDeconParameters()
    {
        Assert.That(_viewModel.Parameters, Is.InstanceOf<MultipleDeconParameters>());
    }
}
