using System;
using System.Collections.ObjectModel;
using System.Linq;
using EngineLayer;
using MassSpectrometry;

namespace GuiFunctions;

public class MultipleDeconParamsViewModel : DeconParamsViewModel
{
    private MultipleDeconParameters _parameters;
    private readonly ObservableCollection<DeconParamsViewModel> _subParameters;
    private DeconvolutionType _selectedAddType = DeconvolutionType.ClassicDeconvolution;

    public ObservableCollection<DeconParamsViewModel> SubParameters => _subParameters;

    public ObservableCollection<DeconvolutionType> AvailableSubTypes { get; } = new()
    {
        DeconvolutionType.ClassicDeconvolution,
        DeconvolutionType.IsoDecDeconvolution
    };

    public DeconvolutionType SelectedAddType
    {
        get => _selectedAddType;
        set
        {
            _selectedAddType = value;
            OnPropertyChanged(nameof(SelectedAddType));
        }
    }

    public override DeconvolutionParameters Parameters
    {
        get => _parameters;
        protected set
        {
            _parameters = (MultipleDeconParameters)value;
            OnPropertyChanged(nameof(Parameters));
        }
    }

    public override string ToString() => "Multiple";

    public MultipleDeconParamsViewModel(MultipleDeconParameters parameters)
    {
        _parameters = parameters;
        _subParameters = new ObservableCollection<DeconParamsViewModel>(
            parameters.Parameters.Select(p => p.ToViewModel()));
    }

    public int SharedMinAssumedChargeState
    {
        get => MinAssumedChargeState;
        set
        {
            MinAssumedChargeState = value;
            foreach (var sub in _subParameters)
                sub.MinAssumedChargeState = value;
        }
    }

    public int SharedMaxAssumedChargeState
    {
        get => MaxAssumedChargeState;
        set
        {
            MaxAssumedChargeState = value;
            foreach (var sub in _subParameters)
                sub.MaxAssumedChargeState = value;
        }
    }

    public new Polarity Polarity
    {
        get => base.Polarity;
        set
        {
            base.Polarity = value;
            foreach (var sub in _subParameters)
                sub.Polarity = value;
        }
    }

    public void AddSubType(DeconvolutionType type, bool isPrecursor = true)
    {
        var defaultParams = type.GetDefaultDeconParams(GlobalVariables.AnalyteType, isPrecursor);
        var subVm = defaultParams.ToViewModel();
        subVm.Polarity = _parameters.Polarity;
        subVm.MinAssumedChargeState = _parameters.MinAssumedChargeState;
        subVm.MaxAssumedChargeState = _parameters.MaxAssumedChargeState;
        _subParameters.Add(subVm);
        RebuildParameters();
    }

    public void RemoveSubType(DeconParamsViewModel vm)
    {
        if (_subParameters.Count <= 1)
            return;

        _subParameters.Remove(vm);
        RebuildParameters();
    }

    private void RebuildParameters()
    {
        _parameters = new MultipleDeconParameters(
            _subParameters.Select(s => s.Parameters),
            _parameters.MinAssumedChargeState,
            _parameters.MaxAssumedChargeState,
            _parameters.Polarity,
            _parameters.AverageResidueModel,
            _parameters.ExpectedIsotopeSpacing);
        OnPropertyChanged(nameof(Parameters));
    }

}
