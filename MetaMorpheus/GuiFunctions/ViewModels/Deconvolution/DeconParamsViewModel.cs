using System;
using System.Diagnostics.CodeAnalysis;
using MassSpectrometry;

namespace GuiFunctions;

/// <summary>
/// Base Deconovolution Parameters view model
/// Used to wrap the DeconvolutionParameters object
/// Contains only shared information between different DeconvolutionParameters
/// </summary>
public abstract class DeconParamsViewModel : BaseViewModel, IEquatable<DeconParamsViewModel>
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


    public override string ToString() => DeconvolutionType.ToString();

    public bool Equals(DeconParamsViewModel other)
    {
        if (ReferenceEquals(null, other)) return false;
        if (ReferenceEquals(this, other)) return true;
        return Equals(Parameters, other.Parameters);
    }

    public override bool Equals(object obj)
    {
        if (ReferenceEquals(null, obj)) return false;
        if (ReferenceEquals(this, obj)) return true;
        if (obj.GetType() != this.GetType()) return false;
        return Equals((DeconParamsViewModel)obj);
    }

    public override int GetHashCode()
    {
        return (Parameters != null ? Parameters.DeconvolutionType.GetHashCode() + Parameters.Polarity.GetHashCode() + MaxAssumedChargeState.GetHashCode() : 0);
    }
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