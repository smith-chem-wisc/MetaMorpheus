using System;
using System.Linq;
using EngineLayer;
using GuiFunctions;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.GuiTests;

[TestFixture]
public class DeconHostViewModelTests
{
    internal ClassicDeconvolutionParameters ClassicPrecursorDeconvolutionParameters = new ClassicDeconvolutionParameters(1, 12, 4, 3);
    internal ClassicDeconvolutionParameters ClassicProductDeconvolutionParameters = new ClassicDeconvolutionParameters(1, 10, 4, 3);

    [Test]
    public void Constructor_DefaultParameters_ShouldInitializeCorrectly()
    {
        // Arrange
        var initialPrecursorParameters = ClassicPrecursorDeconvolutionParameters;
        var initialProductParameters = ClassicProductDeconvolutionParameters;

        // Act
        var viewModel = new DeconHostViewModel(initialPrecursorParameters, initialProductParameters);

        // Assert
        Assert.That(viewModel.UseProvidedPrecursors, Is.False);
        Assert.That(viewModel.DoPrecursorDeconvolution, Is.True);
        Assert.That(viewModel.PrecursorDeconvolutionParametersList, Is.Not.Null);
        Assert.That(viewModel.ProductDeconvolutionParametersList, Is.Not.Null);
        Assert.That(viewModel.PrecursorDeconvolutionParametersList.Any(), Is.True);
        Assert.That(viewModel.ProductDeconvolutionParametersList.Any(), Is.True);
    }

    [Test]
    public void Constructor_WithProvidedParameters_ShouldSetCorrectly()
    {
        // Arrange
        var initialPrecursorParameters = ClassicPrecursorDeconvolutionParameters;
        var initialProductParameters = ClassicProductDeconvolutionParameters;

        // Act
        var viewModel = new DeconHostViewModel(initialPrecursorParameters, initialProductParameters, true, false);

        // Assert
        Assert.That(viewModel.UseProvidedPrecursors, Is.True);
        Assert.That(viewModel.DoPrecursorDeconvolution, Is.False);
        Assert.That(viewModel.PrecursorDeconvolutionParameters, Is.EqualTo(initialPrecursorParameters.ToViewModel()));
        Assert.That(viewModel.ProductDeconvolutionParameters, Is.EqualTo(initialProductParameters.ToViewModel()));
    }

    [Test]
    public void TestDeconHostViewModel_DefaultParameters()
    {
        // Arrange
        DeconHostViewModel viewModel = new DeconHostViewModel(null, null);

        // Act
        var precursorParams = viewModel.PrecursorDeconvolutionParameters;
        var productParams = viewModel.ProductDeconvolutionParameters;

        // Assert
        Assert.That(precursorParams, Is.Not.Null);
        Assert.That(productParams, Is.Not.Null);
        Assert.That(precursorParams.DeconvolutionType, Is.EqualTo(DeconvolutionType.ClassicDeconvolution));
        Assert.That(productParams.DeconvolutionType, Is.EqualTo(DeconvolutionType.ClassicDeconvolution));
    }

    [Test]
    public void PrecursorDeconvolutionParameters_Setter_ShouldTriggerPropertyChanged()
    {
        // Arrange
        var viewModel = new DeconHostViewModel();
        var newParameters = new ClassicDeconParamsViewModel(ClassicPrecursorDeconvolutionParameters);
        bool propertyChangedTriggered = false;
        viewModel.PropertyChanged += (sender, args) =>
        {
            if (args.PropertyName == nameof(viewModel.PrecursorDeconvolutionParameters))
                propertyChangedTriggered = true;
        };

        // Act
        viewModel.PrecursorDeconvolutionParameters = newParameters;

        // Assert
        Assert.That(propertyChangedTriggered, Is.True);
    }

    [Test]
    public void ProductDeconvolutionParameters_Setter_ShouldTriggerPropertyChanged()
    {
        // Arrange
        var viewModel = new DeconHostViewModel();
        var newParameters = new ClassicDeconParamsViewModel(ClassicProductDeconvolutionParameters);
        bool propertyChangedTriggered = false;
        viewModel.PropertyChanged += (sender, args) =>
        {
            if (args.PropertyName == nameof(viewModel.ProductDeconvolutionParameters))
                propertyChangedTriggered = true;
        };

        // Act
        viewModel.ProductDeconvolutionParameters = newParameters;

        // Assert
        Assert.That(propertyChangedTriggered, Is.True);
    }

        
    [Test]
    [NonParallelizable]
    public void TestDeconHostViewModel_GlobalVariables_Proteoform()
    {
        // Arrange
        GlobalVariables.AnalyteType = AnalyteType.Proteoform;
        DeconHostViewModel viewModel = new DeconHostViewModel(null, null);

        // Act
        var precursorParams = viewModel.PrecursorDeconvolutionParameters;
        var productParams = viewModel.ProductDeconvolutionParameters;

        // Assert
        Assert.That(precursorParams, Is.Not.Null);
        Assert.That(productParams, Is.Not.Null);
        Assert.That(precursorParams.DeconvolutionType, Is.EqualTo(DeconvolutionType.ClassicDeconvolution));
        Assert.That(productParams.DeconvolutionType, Is.EqualTo(DeconvolutionType.ClassicDeconvolution));
        Assert.That(precursorParams.Parameters, Is.InstanceOf<ClassicDeconvolutionParameters>());
        Assert.That(productParams.Parameters, Is.InstanceOf<ClassicDeconvolutionParameters>());
        Assert.That(((ClassicDeconvolutionParameters)precursorParams.Parameters).MaxAssumedChargeState, Is.EqualTo(60));
        Assert.That(((ClassicDeconvolutionParameters)productParams.Parameters).MaxAssumedChargeState, Is.EqualTo(10));

        // Revert back to default
        GlobalVariables.AnalyteType = AnalyteType.Peptide;
    }

    [Test]
    [NonParallelizable]
    public void TestDeconHostViewModel_GlobalVariables_Oligo()
    {
        // Arrange
        GlobalVariables.AnalyteType = AnalyteType.Oligo;
        DeconHostViewModel viewModel = new DeconHostViewModel(null, null);

        // Act
        var precursorParams = viewModel.PrecursorDeconvolutionParameters;
        var productParams = viewModel.ProductDeconvolutionParameters;

        // Assert
        Assert.That(precursorParams, Is.Not.Null);
        Assert.That(productParams, Is.Not.Null);
        Assert.That(precursorParams.DeconvolutionType, Is.EqualTo(DeconvolutionType.ClassicDeconvolution));
        Assert.That(productParams.DeconvolutionType, Is.EqualTo(DeconvolutionType.ClassicDeconvolution));
        Assert.That(precursorParams.Parameters, Is.InstanceOf<ClassicDeconvolutionParameters>());
        Assert.That(productParams.Parameters, Is.InstanceOf<ClassicDeconvolutionParameters>());
        Assert.That(((ClassicDeconvolutionParameters)precursorParams.Parameters).MinAssumedChargeState, Is.EqualTo(-20));
        Assert.That(((ClassicDeconvolutionParameters)productParams.Parameters).MinAssumedChargeState, Is.EqualTo(-10));

        // Revert back to default
        GlobalVariables.AnalyteType = AnalyteType.Peptide;
    }

    [Test]
    [NonParallelizable]
    public void TestDeconHostViewModel_GlobalVariables_Unknown()
    {
        // Arrange
        GlobalVariables.AnalyteType = (AnalyteType)(-1);

        // Act & Assert
        Assert.Throws<ArgumentOutOfRangeException>(() =>
        {
            var deconHostViewModel = new DeconHostViewModel(null, null);
        });

        Assert.Throws<ArgumentOutOfRangeException>(() =>
        {
            var deconHostViewModel = new DeconHostViewModel(ClassicPrecursorDeconvolutionParameters, null);
        });

        // Revert back to default
        GlobalVariables.AnalyteType = AnalyteType.Peptide;
    }

    [Test]
    public void TestDisplayDeconSelectionComboBox_MultipleOptions_Precursor()
    {
        // Arrange
        var viewModel = new DeconHostViewModel();
        viewModel.PrecursorDeconvolutionParametersList.Add(ClassicPrecursorDeconvolutionParameters.ToViewModel());
        
        // Act
        var result = viewModel.DisplayDeconSelectionComboBox;

        // Assert
        Assert.That(result, Is.True);
    }

    [Test]
    public void TestDisplayDeconSelectionComboBox_MultipleOptions_Product()
    {
        // Arrange
        var viewModel = new DeconHostViewModel();
        viewModel.ProductDeconvolutionParametersList.Add(ClassicProductDeconvolutionParameters.ToViewModel());

        // Act
        var result = viewModel.DisplayDeconSelectionComboBox;

        // Assert
        Assert.That(result, Is.True);
    }

    [Test]
    public void TestDisplayDeconSelectionComboBox_SingleOption()
    {
        // Arrange
        var viewModel = new DeconHostViewModel();
       
        // Act
        var result = viewModel.DisplayDeconSelectionComboBox;

        // Assert
        Assert.That(result, Is.False);
    }
}