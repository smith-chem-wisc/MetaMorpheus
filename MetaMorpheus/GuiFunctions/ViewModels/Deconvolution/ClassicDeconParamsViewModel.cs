using MassSpectrometry;

namespace GuiFunctions;

public class ClassicDeconParamsViewModel : DeconParamsViewModel
{
    private ClassicDeconvolutionParameters _parameters;
    public override DeconvolutionParameters Parameters
    {
        get => _parameters;
        protected set
        {
            _parameters = (ClassicDeconvolutionParameters)value;
            OnPropertyChanged(nameof(Parameters));
        }
    }

    public ClassicDeconParamsViewModel(ClassicDeconvolutionParameters parameters)
    {
        _parameters = parameters;
    }

    public double DeconvolutionTolerancePpm
    {
        get => _parameters.DeconvolutionTolerancePpm;
        set
        {
            _parameters.DeconvolutionTolerancePpm = value;
            OnPropertyChanged(nameof(DeconvolutionTolerancePpm));
        }
    }

    public double IntensityRatioLimit
    {
        get => _parameters.IntensityRatioLimit;
        set
        {
            _parameters.IntensityRatioLimit = value;
            OnPropertyChanged(nameof(IntensityRatioLimit));
        }
    }
}