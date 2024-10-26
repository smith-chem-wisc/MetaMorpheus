#nullable enable
using System;
using System.Collections.ObjectModel;
using System.Diagnostics.CodeAnalysis;
using EngineLayer;
using MassSpectrometry;

namespace GuiFunctions;

public class DeconHostViewModel : BaseViewModel
{

    public DeconHostViewModel(DeconvolutionParameters? initialPrecursorParameters = null, DeconvolutionParameters? initialFragmentParameters = null,
        bool useProvidedPrecursor = false, bool deconvolutePrecursors = true)
    {
        // Always the same
        DeconvolutionTypes = new ObservableCollection<DeconvolutionType>
        {
            DeconvolutionType.ClassicDeconvolution
        };

        // Parameter Dependent
        UseProvidedPrecursors = useProvidedPrecursor;
        DoPrecursorDeconvolution = deconvolutePrecursors;

        initialPrecursorParameters ??= GlobalVariables.AnalyteType switch
        {
            "Peptide" => new ClassicDeconvolutionParameters(1, 12, 4, 3),
            "Proteoform" => new ClassicDeconvolutionParameters(1, 60, 4, 3),
            "Oligo" => new ClassicDeconvolutionParameters(-20, -1, 4, 3),
            _ => throw new ArgumentOutOfRangeException()
        };
        PrecursorDeconvolutionParameters = initialPrecursorParameters.ToViewModel();

        initialFragmentParameters ??= GlobalVariables.AnalyteType switch
        {
            "Peptide" => new ClassicDeconvolutionParameters(1, 12, 4, 3),
            "Proteoform" => new ClassicDeconvolutionParameters(1, 60, 4, 3),
            "Oligo" => new ClassicDeconvolutionParameters(-20, -1, 4, 3),
            _ => throw new ArgumentOutOfRangeException()
        };
        ProductDeconvolutionParameters = initialFragmentParameters.ToViewModel();
    }

    #region Common Parameters

    private bool _useProvidedPrecursors;
    public bool UseProvidedPrecursors
    {
        get => _useProvidedPrecursors;
        set
        {
            _useProvidedPrecursors = value;
            OnPropertyChanged(nameof(UseProvidedPrecursors));
        }
    }

    private bool _doPrecursorDeconvolution;
    public bool DoPrecursorDeconvolution
    {
        get => _doPrecursorDeconvolution;
        set
        {
            _doPrecursorDeconvolution = value;
            OnPropertyChanged(nameof(DoPrecursorDeconvolution));
        }
    }

    #endregion

    public ObservableCollection<DeconvolutionType> DeconvolutionTypes { get; private set; }

    private DeconvolutionType _selectedDeconvolutionType;
    public DeconvolutionType SelectedDeconvolutionType
    {
        get => _selectedDeconvolutionType;
        set
        {
            _selectedDeconvolutionType = value;
            OnPropertyChanged(nameof(SelectedDeconvolutionType));
        }
    }


    private DeconParamsViewModel _precursorDeconvolutionParameters;
    public DeconParamsViewModel PrecursorDeconvolutionParameters
    {
        get => _precursorDeconvolutionParameters;
        set
        {
            _precursorDeconvolutionParameters = value;
            OnPropertyChanged(nameof(PrecursorDeconvolutionParameters));
        }
    }

    private DeconParamsViewModel _productDeconvolutionParameters;
    public DeconParamsViewModel ProductDeconvolutionParameters
    {
        get => _productDeconvolutionParameters;
        set
        {
            _productDeconvolutionParameters = value;
            OnPropertyChanged(nameof(ProductDeconvolutionParameters));
        }
    }
}



[ExcludeFromCodeCoverage] // Model used only for visualizing the view
public class DeconHostModel : DeconHostViewModel
{
    public static DeconHostModel Instance => new DeconHostModel();

    public DeconHostModel() : base (DeconParamsModel.Instance.Parameters, DeconParamsModel.Instance.Parameters)
    {
        UseProvidedPrecursors = false;
        DoPrecursorDeconvolution = true;
    }
}