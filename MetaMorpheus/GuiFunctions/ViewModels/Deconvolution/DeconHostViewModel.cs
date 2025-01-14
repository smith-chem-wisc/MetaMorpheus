#nullable enable
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Easy.Common.Extensions;
using EngineLayer;
using MassSpectrometry;

namespace GuiFunctions;

/// <summary>
/// This class holds all of the information in the Deconvolution tab of the GUI
/// One instance will be created per Task Window
///
/// The Task window will populate this view model with the appropriate parameters from <see cref="CommonParameters"/>
/// The user can then modify these parameters as needed via the gui
/// The Task window will use the PrecursorDeconvolutionParameters and teh ProductDeconvolutionParameters to create a new <see cref="CommonParameters"/> object
/// </summary>
public class DeconHostViewModel : BaseViewModel
{
    /// <summary>
    /// This is where default deconvolution parameters are set for GUI display
    /// </summary>
    /// <param name="initialPrecursorParameters">precursor params to display first</param>
    /// <param name="initialProductParameters">product params to display first</param>
    /// <param name="useProvidedPrecursor"></param>
    /// <param name="deconvolutePrecursors"></param>
    /// <exception cref="ArgumentOutOfRangeException"></exception>
    public DeconHostViewModel(DeconvolutionParameters? initialPrecursorParameters = null, DeconvolutionParameters? initialProductParameters = null,
        bool useProvidedPrecursor = false, bool deconvolutePrecursors = true)
    {
        UseProvidedPrecursors = useProvidedPrecursor;
        DoPrecursorDeconvolution = deconvolutePrecursors;

        // Order matters here, construct the lists before setting the selected parameters
        PrecursorDeconvolutionParametersList = new ObservableCollection<DeconParamsViewModel>();
        ProductDeconvolutionParametersList = new ObservableCollection<DeconParamsViewModel>();

        // populate the lists by adding the default parameters for each deconvolution type or the provided parameters
        foreach (var deconType in Enum.GetValues<DeconvolutionType>())
        {
            switch (deconType)
            {
                case DeconvolutionType.ExampleNewDeconvolutionTemplate:
                    continue;

                case DeconvolutionType.ClassicDeconvolution:

                    // Precursor
                    if (initialPrecursorParameters is { DeconvolutionType: DeconvolutionType.ClassicDeconvolution })
                        PrecursorDeconvolutionParametersList.Add(initialPrecursorParameters.ToViewModel());
                    else
                    {
                        var toAdd = GlobalVariables.AnalyteType switch
                        {
                            AnalyteType.Peptide => new ClassicDeconvolutionParameters(1, 12, 4, 3),
                            AnalyteType.Proteoform => new ClassicDeconvolutionParameters(1, 60, 4, 3),
                            AnalyteType.Oligo => new ClassicDeconvolutionParameters(-20, -1, 4, 3),
                            _ => throw new ArgumentOutOfRangeException()
                        };
                        PrecursorDeconvolutionParametersList.Add(toAdd.ToViewModel());
                    }
                    
                    // Product
                    if (initialProductParameters is { DeconvolutionType: DeconvolutionType.ClassicDeconvolution })
                        ProductDeconvolutionParametersList.Add(initialProductParameters.ToViewModel());
                    else
                    {
                        var toAdd = GlobalVariables.AnalyteType switch
                        {
                            AnalyteType.Peptide => new ClassicDeconvolutionParameters(1, 10, 4, 3),
                            AnalyteType.Proteoform => new ClassicDeconvolutionParameters(1, 10, 4, 3),
                            AnalyteType.Oligo => new ClassicDeconvolutionParameters(-10, -1, 4, 3),
                            _ => throw new ArgumentOutOfRangeException()
                        };
                        ProductDeconvolutionParametersList.Add(toAdd.ToViewModel());
                    }

                    break;

                case DeconvolutionType.IsoDecDeconvolution:

                    // Precursor
                    if (initialPrecursorParameters is { DeconvolutionType: DeconvolutionType.IsoDecDeconvolution })
                        PrecursorDeconvolutionParametersList.Add(initialPrecursorParameters.ToViewModel());
                    else
                    {
                        var toAdd = GlobalVariables.AnalyteType switch
                        {
                            AnalyteType.Peptide => new IsoDecDeconvolutionParameters(),
                            AnalyteType.Proteoform => new IsoDecDeconvolutionParameters(),
                            AnalyteType.Oligo => new IsoDecDeconvolutionParameters(Polarity.Negative),
                            _ => throw new ArgumentOutOfRangeException()
                        };
                        PrecursorDeconvolutionParametersList.Add(toAdd.ToViewModel());
                    }

                    // Product
                    if (initialProductParameters is { DeconvolutionType: DeconvolutionType.IsoDecDeconvolution })
                        ProductDeconvolutionParametersList.Add(initialProductParameters.ToViewModel());
                    else
                    {
                        var toAdd = GlobalVariables.AnalyteType switch
                        {
                            AnalyteType.Peptide => new IsoDecDeconvolutionParameters(reportMultipleMonoisos: false),
                            AnalyteType.Proteoform => new IsoDecDeconvolutionParameters(reportMultipleMonoisos: false),
                            AnalyteType.Oligo => new IsoDecDeconvolutionParameters(Polarity.Negative, reportMultipleMonoisos: false),
                            _ => throw new ArgumentOutOfRangeException()
                        };
                        ProductDeconvolutionParametersList.Add(toAdd.ToViewModel());
                    }

                    break;

                default: // This will only be hit if a new deconvolution type is added to mzlib and not handled here
                    throw new ArgumentOutOfRangeException();
            }
        }

        // If deconvolution parameters are not set, default to MetaMorpheus defaults
        PrecursorDeconvolutionParameters = initialPrecursorParameters is null 
            ? PrecursorDeconvolutionParametersList.First(x => x.DeconvolutionType == DeconvolutionType.ClassicDeconvolution) 
            : PrecursorDeconvolutionParametersList.First(x => x.Parameters == initialPrecursorParameters);

        ProductDeconvolutionParameters = initialProductParameters is null 
            ? ProductDeconvolutionParametersList.First(x => x.DeconvolutionType == DeconvolutionType.ClassicDeconvolution)
            : ProductDeconvolutionParametersList.First(x => x.Parameters == initialProductParameters);
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

    public void SetAllPrecursorMaxChargeState(int newMaxCharge)
    {
        foreach (var precursorParams in PrecursorDeconvolutionParametersList)
        {
            precursorParams.MaxAssumedChargeState = newMaxCharge;
        }
        OnPropertyChanged(nameof(PrecursorDeconvolutionParametersList));
    }

    public void SetAllProductMaxChargeState(int newMaxCharge)
    {
        foreach (var productParams in ProductDeconvolutionParametersList)
        {
            productParams.MaxAssumedChargeState = newMaxCharge;
        }
        OnPropertyChanged(nameof(ProductDeconvolutionParametersList));
    }

    #endregion

    /// <summary>
    /// All of the possible precursor deconvolution parameters that can be selected
    /// 
    /// Their ToString() method sets the name of the combo box
    /// Stored in memory while task window is open, only the selected one is used for <See cref="CommonParameters"/>
    /// This enables the user to set parameters, switch to another, and switch back without losing their settings
    /// </summary>
    public ObservableCollection<DeconParamsViewModel> PrecursorDeconvolutionParametersList { get; protected set; }
    private DeconParamsViewModel? _precursorDeconvolutionParameters;

    /// <summary>
    /// The selected precursor deconvolution parameters
    /// </summary>
    public DeconParamsViewModel PrecursorDeconvolutionParameters
    {
        get => _precursorDeconvolutionParameters!;
        set
        {
            _precursorDeconvolutionParameters = value;
            OnPropertyChanged(nameof(PrecursorDeconvolutionParameters));
        }
    }

    /// <summary>
    /// All of the possible product deconvolution parameters that can be selected
    /// 
    /// Their ToString() method sets the name of the combo box
    /// Stored in memory while task window is open, only the selected one is used for <See cref="CommonParameters"/>
    /// This enables the user to set parameters, switch to another, and switch back without losing their settings
    /// </summary>
    public ObservableCollection<DeconParamsViewModel> ProductDeconvolutionParametersList { get; protected set; }
    private DeconParamsViewModel? _productDeconvolutionParameters;

    /// <summary>
    /// The selected product deconvolution parameters
    /// </summary>
    public DeconParamsViewModel ProductDeconvolutionParameters
    {
        get => _productDeconvolutionParameters!;
        set
        {
            _productDeconvolutionParameters = value;
            OnPropertyChanged(nameof(ProductDeconvolutionParameters));
        }
    }

    /// <summary>
    /// Hides the decon type selection combo box if only one options is present
    /// </summary>
    public bool DisplayDeconSelectionComboBox => PrecursorDeconvolutionParametersList.Count > 1 || ProductDeconvolutionParametersList.Count > 1;
}

[ExcludeFromCodeCoverage] // Model used only for visualizing the view in visual studio
public class DeconHostModel : DeconHostViewModel
{
    public static DeconHostModel Instance => new();

    public DeconHostModel() : base (DeconParamsModel.Instance.Parameters, DeconParamsModel.Instance.Parameters)
    {
        UseProvidedPrecursors = false;
        DoPrecursorDeconvolution = true;
    }
}