using GuiFunctions;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.GuiTests;

[TestFixture]
public class DeconParamsViewModelTest
{
    private class TestDeconParamsViewModel : DeconParamsViewModel
    {
        public sealed override DeconvolutionParameters Parameters { get; protected set; }

        public TestDeconParamsViewModel(DeconvolutionParameters parameters)
        {
            Parameters = parameters;
        }
    }

    private DeconvolutionParameters _negativeModeParameters;
    private TestDeconParamsViewModel _negativeModeViewModel;
    private DeconvolutionParameters _positiveModeParameters; 
    private TestDeconParamsViewModel _positiveModeViewModel;

    [SetUp]
    public void SetUp()
    {
        _negativeModeParameters = new ClassicDeconvolutionParameters(-20, -2, 5, 3, Polarity.Negative);
        _negativeModeViewModel = new TestDeconParamsViewModel(_negativeModeParameters);

        _positiveModeParameters = new ClassicDeconvolutionParameters(2, 20, 5, 3, Polarity.Positive);
        _positiveModeViewModel = new TestDeconParamsViewModel(_positiveModeParameters);
    }

    [TestCase(0, false, TestName = "MinAssumedChargeState_SetToZero_DoesNotChange")]
    [TestCase(-1, false, TestName = "MinAssumedChargeState_SetToNegativeInPositiveMode_DoesNotChange")]
    [TestCase(21, false, TestName = "MinAssumedChargeState_SetAboveMaxInPositiveMode_DoesNotChange")]
    [TestCase(3, true, TestName = "MinAssumedChargeState_SetToValidPositiveValue_Changes")]
    public void TestMinAssumedChargeState_PositiveMode(int newValue, bool shouldChange)
    {
        _positiveModeViewModel.MinAssumedChargeState = newValue;
        if (shouldChange)
        {
            Assert.That(_positiveModeViewModel.MinAssumedChargeState, Is.EqualTo(newValue));
        }
        else
        {
            Assert.That(_positiveModeViewModel.MinAssumedChargeState, Is.EqualTo(2));
        }
    }

    [TestCase(0, false, TestName = "MaxAssumedChargeState_SetToZero_DoesNotChange")]
    [TestCase(-4, false, TestName = "MaxAssumedChargeState_SetToNegativeInPositiveMode_DoesNotChange")]
    [TestCase(1, false, TestName = "MaxAssumedChargeState_SetBelowMinInPositiveMode_DoesNotChange")]
    [TestCase(25, true, TestName = "MaxAssumedChargeState_SetToValidPositiveValue_Changes")]
    public void TestMaxAssumedChargeState_PositiveMode(int newValue, bool shouldChange)
    {
        _positiveModeViewModel.MaxAssumedChargeState = newValue;
        if (shouldChange)
        {
            Assert.That(_positiveModeViewModel.MaxAssumedChargeState, Is.EqualTo(newValue));
        }
        else
        {
            Assert.That(_positiveModeViewModel.MaxAssumedChargeState, Is.EqualTo(20));
        }
    }

    [TestCase(0, false, TestName = "MinAssumedChargeState_SetToZero_DoesNotChange")]
    [TestCase(1, false, TestName = "MinAssumedChargeState_SetToPositiveInNegativeMode_DoesNotChange")]
    [TestCase(-1, false, TestName = "MinAssumedChargeState_SetAboveMaxInNegativeMode_DoesNotChange")]
    [TestCase(-21, true, TestName = "MinAssumedChargeState_SetToValidNegativeValue_Changes")]
    public void TestMinAssumedChargeState_NegativeMode(int newValue, bool shouldChange)
    {
        _negativeModeViewModel.MinAssumedChargeState = newValue;
        if (shouldChange)
        {
            Assert.That(_negativeModeViewModel.MinAssumedChargeState, Is.EqualTo(newValue));
        }
        else
        {
            Assert.That(_negativeModeViewModel.MinAssumedChargeState, Is.EqualTo(-20));
        }
    }

    [TestCase(0, false, TestName = "MaxAssumedChargeState_SetToZero_DoesNotChange")]
    [TestCase(1, false, TestName = "MaxAssumedChargeState_SetToPositiveInNegativeMode_DoesNotChange")]
    [TestCase(-21, false, TestName = "MaxAssumedChargeState_SetBelowMinInNegativeMode_DoesNotChange")]
    [TestCase(-1, true, TestName = "MaxAssumedChargeState_SetToValidNegativeValue_Changes")]
    public void TestMaxAssumedChargeState_NegativeMode(int newValue, bool shouldChange)
    {
        _negativeModeViewModel.MaxAssumedChargeState = newValue;
        if (shouldChange)
        {
            Assert.That(_negativeModeViewModel.MaxAssumedChargeState, Is.EqualTo(newValue));
        }
        else
        {
            Assert.That(_negativeModeViewModel.MaxAssumedChargeState, Is.EqualTo(-2));
        }
    }

    [TestCase(Polarity.Positive, Polarity.Positive, 5, 12, TestName = "Polarity_SwitchFromPositiveToPositive_NoChange")]
    [TestCase(Polarity.Negative, Polarity.Negative, -25, -5, TestName = "Polarity_SwitchFromNegativeToNegative_NoChange")]
    [TestCase(Polarity.Positive, Polarity.Negative, -20, -1, TestName = "Polarity_SwitchFromPositiveToNegative_UpdatesChargeStates")]
    [TestCase(Polarity.Negative, Polarity.Positive, 1, 12, TestName = "Polarity_SwitchFromNegativeToPositive_UpdatesChargeStates")]
    public void TestPolaritySwitch(Polarity initialPolarity, Polarity newPolarity, int expectedMinChargeState, int expectedMaxChargeState)
    {
        var parameters = new ClassicDeconvolutionParameters(initialPolarity == Polarity.Positive ? 5 : -25, initialPolarity == Polarity.Positive ? 12 : -5, 5, 3, initialPolarity);
        var viewModel = new TestDeconParamsViewModel(parameters);

        viewModel.Polarity = newPolarity;

        Assert.That(viewModel.Polarity, Is.EqualTo(newPolarity));
        Assert.That(viewModel.MinAssumedChargeState, Is.EqualTo(expectedMinChargeState));
        Assert.That(viewModel.MaxAssumedChargeState, Is.EqualTo(expectedMaxChargeState));
    }

    [Test]
    public void TestPolaritySwitchWithPreviousValues_PositiveMode()
    {
        // Initial setup with positive polarity
        var parameters = new ClassicDeconvolutionParameters(1, 12, 5, 3, Polarity.Positive);
        var viewModel = new TestDeconParamsViewModel(parameters);

        // Set previous values
        viewModel.MinAssumedChargeState = 5;
        viewModel.MaxAssumedChargeState = 10;

        // Change to negative polarity
        viewModel.Polarity = Polarity.Negative;

        // Check that defaults got set properly
        Assert.That(viewModel.MinAssumedChargeState, Is.EqualTo(-20));
        Assert.That(viewModel.MaxAssumedChargeState, Is.EqualTo(-1));

        // Change back to positive polarity
        viewModel.Polarity = Polarity.Positive;

        // Assert that previous values are restored correctly
        Assert.That(viewModel.MinAssumedChargeState, Is.EqualTo(5));
        Assert.That(viewModel.MaxAssumedChargeState, Is.EqualTo(10));
    }

    [Test]
    public void TestPolaritySwitchWithPreviousValues_NegativeMode()
    {
        // Initial setup with negative polarity
        var parameters = new ClassicDeconvolutionParameters(-20, -1, 5, 3, Polarity.Negative);
        var viewModel = new TestDeconParamsViewModel(parameters);

        // Set previous values
        viewModel.MinAssumedChargeState = -10;
        viewModel.MaxAssumedChargeState = -5;

        // Change to positive polarity
        viewModel.Polarity = Polarity.Positive;

        // Check that defaults got set properly
        Assert.That(viewModel.MinAssumedChargeState, Is.EqualTo(1));
        Assert.That(viewModel.MaxAssumedChargeState, Is.EqualTo(12));

        // Change back to negative polarity
        viewModel.Polarity = Polarity.Negative;

        // Assert that previous values are restored correctly
        Assert.That(viewModel.MinAssumedChargeState, Is.EqualTo(-10));
        Assert.That(viewModel.MaxAssumedChargeState, Is.EqualTo(-5));
    }

    [Test]
    public void TestToString()
    {
        var parameters = new ClassicDeconvolutionParameters(1, 12, 5, 3, Polarity.Positive);
        var viewModel = new TestDeconParamsViewModel(parameters);
        Assert.That(viewModel.ToString(), Is.EqualTo(parameters.DeconvolutionType.ToString()));
    }

    [Test]
    public void TestEquals_SameObject()
    {
        var parameters = new ClassicDeconvolutionParameters(1, 12, 5, 3, Polarity.Positive);
        var viewModel = new TestDeconParamsViewModel(parameters);
        Assert.That(viewModel.Equals(viewModel), Is.True);
    }

    [Test]
    public void TestEquals_NullObject()
    {
        var parameters = new ClassicDeconvolutionParameters(1, 12, 5, 3, Polarity.Positive);
        var viewModel = new TestDeconParamsViewModel(parameters);
        Assert.That(viewModel.Equals(null), Is.False);
    }

    [Test]
    public void TestEquals_DifferentType()
    {
        var parameters = new ClassicDeconvolutionParameters(1, 12, 5, 3, Polarity.Positive);
        var viewModel = new TestDeconParamsViewModel(parameters);
        Assert.That(viewModel.Equals(new object()), Is.False);
    }

    [Test]
    [TestCase(1, 15, Polarity.Negative, TestName = "MinChargeDiffers_NotEqual")]
    [TestCase(2, 12, Polarity.Negative, TestName = "MaxChargeDiffers_NotEqual")]
    [TestCase(2, 15, Polarity.Positive, TestName = "PolarityDiffers_NotEqual")]
    public void TestEquals_DifferentParameters(int minCharge, int maxCharge, Polarity polarity)
    {
        var parameters1 = new ClassicDeconvolutionParameters(minCharge, maxCharge, 5, 3, polarity);
        var viewModel1 = new TestDeconParamsViewModel(parameters1);

        var parameters2 = new ClassicDeconvolutionParameters(2, 15, 5, 3, Polarity.Negative);
        var viewModel2 = new TestDeconParamsViewModel(parameters2);
        Assert.That(viewModel1.Equals(viewModel2), Is.False);
    }

    [Test]
    public void TestEquals_SameParameters()
    {
        var parameters1 = new ClassicDeconvolutionParameters(1, 12, 5, 3, Polarity.Positive);
        var viewModel1 = new TestDeconParamsViewModel(parameters1);
        var parameters2 = new ClassicDeconvolutionParameters(1, 12, 5, 3, Polarity.Positive);
        var viewModel2 = new TestDeconParamsViewModel(parameters2);
        Assert.That(viewModel1.Equals(viewModel2), Is.True);
    }

    [Test]
    public void TestGetHashCode()
    {
        var parameters = new ClassicDeconvolutionParameters(1, 12, 5, 3, Polarity.Positive);
        var viewModel = new TestDeconParamsViewModel(parameters);
        var expectedHashCode = parameters.DeconvolutionType.GetHashCode() + parameters.Polarity.GetHashCode() + parameters.MaxAssumedChargeState.GetHashCode();
        Assert.That(viewModel.GetHashCode(), Is.EqualTo(expectedHashCode));
    }
}