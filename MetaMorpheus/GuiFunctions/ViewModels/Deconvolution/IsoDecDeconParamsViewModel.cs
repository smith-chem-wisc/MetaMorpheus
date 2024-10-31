using MassSpectrometry;

namespace GuiFunctions;

public sealed class IsoDecDeconParamsViewModel : DeconParamsViewModel
{
    private IsoDecDeconvolutionParameters _parameters;
    public override DeconvolutionParameters Parameters
    {
        get => _parameters;
        protected set
        {
            _parameters = (IsoDecDeconvolutionParameters)value;
            OnPropertyChanged(nameof(Parameters));
        }
    }

    public IsoDecDeconParamsViewModel(IsoDecDeconvolutionParameters parameters)
    {
        Parameters = parameters;
    }

    public int PhaseRes
    {
        get => _parameters.PhaseRes;
        set
        {
            if (value != 4 && value != 8)
                return;

            _parameters.PhaseRes = value;
            OnPropertyChanged(nameof(PhaseRes));
        }
    }

    public float CssThreshold
    {
        get => _parameters.CssThreshold;
        set
        {
            if (value is > 1 or < 0)
                return;

            _parameters.CssThreshold = value;
            OnPropertyChanged(nameof(CssThreshold));
        }
    }

    public float MatchTolerance
    {
        get => _parameters.MatchTolerance;
        set
        {
            if (value <= 0)
                return;

            _parameters.MatchTolerance = value;
            OnPropertyChanged(nameof(MatchTolerance));
        }
    }

    public int MaximumShift
    {
        get => _parameters.MaxShift;
        set
        {
            if (value < 0)
                return;

            _parameters.MaxShift = value;
            OnPropertyChanged(nameof(MaximumShift));
        }
    }

    public float MzWindowForIsotopeDistributionMinimum
    {
        get => _parameters.MzWindow[0];
        set
        {
            if (value > 0)
                return;

            _parameters.MzWindow[0] = value;
            OnPropertyChanged(nameof(MzWindowForIsotopeDistributionMinimum));
        }
    } 

    public float MzWindowForIsotopeDistributionMaximum
    {
        get => _parameters.MzWindow[1];
        set
        {
            if (value < 0)
                return;

            _parameters.MzWindow[1] = value;
            OnPropertyChanged(nameof(MzWindowForIsotopeDistributionMaximum));
        }
    }

    public int KnockdownRounds
    {
        get => _parameters.KnockdownRounds;
        set
        {
            if (value <0)
                return;

            _parameters.KnockdownRounds = value;
            OnPropertyChanged(nameof(KnockdownRounds));
        }
    }

    public float MinAreaCovered
    {
        get => _parameters.MinAreaCovered;
        set
        {
            if (value is > 1 or < 0)
                return;

            _parameters.MinAreaCovered = value;
            OnPropertyChanged(nameof(MinAreaCovered));
        }
    }

    public float DataThreshold
    {
        get => _parameters.DataThreshold;
        set
        {
            if (value is > 1 or < 0)
                return;

            _parameters.DataThreshold = value;
            OnPropertyChanged(nameof(DataThreshold));
        }
    }

    public bool ReportMultipleMonoisotopicMasses
    {
        get => _parameters.ReportMulitpleMonoisos;
        set
        {
            _parameters.ReportMulitpleMonoisos = value;
            OnPropertyChanged(nameof(ReportMultipleMonoisotopicMasses));
        }
    }

    public override string ToString() => "IsoDec";
}