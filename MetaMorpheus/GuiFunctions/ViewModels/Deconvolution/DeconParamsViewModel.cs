using System.Diagnostics.CodeAnalysis;
using MassSpectrometry;

namespace GuiFunctions;

public abstract class DeconParamsViewModel : BaseViewModel
{
    public DeconvolutionType DeconvolutionType => Parameters.DeconvolutionType;
    public abstract DeconvolutionParameters Parameters { get; protected set; }

    #region Wrapped Deconvolution Parameters

    public int MinAssumedChargeState
    {
        get => Parameters.MinAssumedChargeState;
        set
        {
            Parameters.MinAssumedChargeState = value;
            OnPropertyChanged(nameof(MinAssumedChargeState));
        }
    }

    public int MaxAssumedChargeState
    {
        get => Parameters.MaxAssumedChargeState;
        set
        {
            Parameters.MaxAssumedChargeState = value;
            OnPropertyChanged(nameof(MaxAssumedChargeState));
        }
    }

    public Polarity Polarity
    {
        get => Parameters.Polarity;
        set
        {
            Parameters.Polarity = value;
            OnPropertyChanged(nameof(Polarity));
        }
    }

    #endregion

}


[ExcludeFromCodeCoverage] // Model used only for visualizing the view
public class DeconParamsModel : DeconParamsViewModel
{
    public static DeconParamsModel Instance => new DeconParamsModel();

    public sealed override DeconvolutionParameters Parameters { get; protected set; }

    public DeconParamsModel()
    {
        Parameters = new ClassicDeconvolutionParameters(1, 20, 5, 3, Polarity.Negative);
    }
}