using System;
using System.Collections.Generic;
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

    private int _previousMinAssumedChargeState;
    private int _previousMaxAssumedChargeState;

    public DeconvolutionType DeconvolutionType => Parameters.DeconvolutionType;
    public abstract DeconvolutionParameters Parameters { get; protected set; }


    /// <summary>
    /// Gets or sets the minimum assumed charge state.
    /// Ensures the value is within valid range based on the polarity.
    /// </summary>
    public int MinAssumedChargeState
    {
        get => Parameters.MinAssumedChargeState;
        set
        {
            if (value == 0)
                return;
            switch (Polarity)
            {
                case Polarity.Positive when value < 0:
                case Polarity.Positive when value > MaxAssumedChargeState:
                case Polarity.Negative when value > 0:
                case Polarity.Negative when value > MaxAssumedChargeState:
                    return;

                default:
                    _previousMinAssumedChargeState = Parameters.MinAssumedChargeState;
                    Parameters.MinAssumedChargeState = value;
                    OnPropertyChanged(nameof(MinAssumedChargeState));
                    break;
            }
        }
    }


    ///<summary>
    /// Gets or sets the maximum assumed charge state.
    /// Ensures the value is within valid range based on the polarity.
    /// </summary>
    public int MaxAssumedChargeState
    {
        get => Parameters.MaxAssumedChargeState;
        set
        {
            if (value == 0)
                return;
            switch (Polarity)
            {
                case Polarity.Positive when value < 0:
                case Polarity.Positive when value < MinAssumedChargeState:
                case Polarity.Negative when value > 0:
                case Polarity.Negative when value < MinAssumedChargeState:
                    return;

                default:
                    _previousMaxAssumedChargeState = Parameters.MaxAssumedChargeState;
                    Parameters.MaxAssumedChargeState = value;
                    OnPropertyChanged(nameof(MaxAssumedChargeState));
                    break;

            }
        }
    }

    // TODO: When RNA is added to MetaMorpheus, add the polarity to the GUI
    public Polarity Polarity
    {
        get => Parameters.Polarity;
        set
        {
            // Update based upon the new polarity
            Parameters.Polarity = value;
            OnPropertyChanged(nameof(Polarity));

            // Update min and max charge to be valid with the new polarity
            // switch from negative to positive
            if (value == Polarity.Positive)
            {
                if (MaxAssumedChargeState < 0)
                {
                    if (_previousMaxAssumedChargeState > 0)
                        MaxAssumedChargeState = _previousMaxAssumedChargeState;
                    else
                        MaxAssumedChargeState = 12;
                }

                if (MinAssumedChargeState < 1)
                {
                    if (_previousMinAssumedChargeState > 0 && _previousMinAssumedChargeState < MaxAssumedChargeState)
                        MinAssumedChargeState = _previousMinAssumedChargeState;
                    else
                        MinAssumedChargeState = 1;
                }
            }

            // switch from positive to negative
            if (value == Polarity.Negative)
            {
                if (MinAssumedChargeState > 0)
                {
                    if (_previousMinAssumedChargeState < 0)
                        MinAssumedChargeState = _previousMinAssumedChargeState;
                    else
                        MinAssumedChargeState = -20;
                }

                if (MaxAssumedChargeState > 0)
                {
                    if (_previousMaxAssumedChargeState < 0 && _previousMaxAssumedChargeState > MinAssumedChargeState)
                        MaxAssumedChargeState = _previousMaxAssumedChargeState;
                    else
                        MaxAssumedChargeState = -1;
                }
            }

        }
    }
    
    public override string ToString() => DeconvolutionType.ToString();

    public bool Equals(DeconParamsViewModel other)
    {
        if (ReferenceEquals(null, other)) return false;
        if (ReferenceEquals(this, other)) return true;

        if (Parameters.DeconvolutionType != other.Parameters.DeconvolutionType)
            return false;
        if (Parameters.Polarity != other.Parameters.Polarity)
            return false;
        if (Parameters.MinAssumedChargeState != other.Parameters.MinAssumedChargeState)
            return false;
        if (Parameters.MaxAssumedChargeState != other.Parameters.MaxAssumedChargeState)
            return false;
        return true;
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


[ExcludeFromCodeCoverage] // Model used only for visualizing the view in visual studio
public class DeconParamsModel : DeconParamsViewModel
{
    public static DeconParamsModel Instance => new DeconParamsModel();

    public sealed override DeconvolutionParameters Parameters { get; protected set; }

    public DeconParamsModel()
    {
        Parameters = new ClassicDeconvolutionParameters(1, 20, 5, 3, Polarity.Negative);
    }
}