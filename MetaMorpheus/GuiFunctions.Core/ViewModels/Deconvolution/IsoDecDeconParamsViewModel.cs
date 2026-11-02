using MassSpectrometry;

namespace GuiFunctions;

// Todo: IsoDec has the ability to report all plausible MonoIsotopic masses for a peak (missed monos). 
// This behavior is currently suppressed due to the complexities of interacting with our mass difference acceptor. 
//      The Isotpic Envelopes that IsoDec produces are labeled with an ID.
//      This ID can be used in PostSearchAnalysis to ensure that no spectral match can be made from two deconvolution results of the same species.

public sealed class IsoDecDeconParamsViewModel : DeconParamsViewModel
{
    public static IsoDecDeconParamsViewModel Instance => new (new IsoDecDeconvolutionParameters());
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
        // Todo: remove this and reconcile the Missed Monoisotopics
        parameters.ReportMulitpleMonoisos = false;

        Parameters = parameters;
    }

    public bool PhaseResIsFour
    {
        get => _parameters.PhaseRes == 4;
        set
        {
            if (value)
            {
                _parameters.PhaseRes = 4;
                OnPropertyChanged(nameof(PhaseResIsFour));
            }
            else
            {
                _parameters.PhaseRes = 8;
                OnPropertyChanged(nameof(PhaseResIsFour));
            }
        }
    }

    /// <summary>
    /// Minimum cosine similarity score for isotope distribution
    /// </summary>
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

    /// <summary>
    /// Match Tolerance for peak detection in ppm
    /// </summary>
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

    /// <summary>
    /// Maximum shift allowed for isotope distribution
    /// </summary>
    public int MaximumShift
    {
        get => _parameters.MaxShift;
        set
        {
            if (value is < 0 or > 8)
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

    /// <summary>
    /// Number of knockdown rounds
    /// </summary>
    public int KnockdownRounds
    {
        get => _parameters.KnockdownRounds;
        set
        {
            if (value < 0)
                return;

            _parameters.KnockdownRounds = value;
            OnPropertyChanged(nameof(KnockdownRounds));
        }
    }

    /// <summary>
    /// Minimum area covered by isotope distribution. Use in or with css_thresh
    /// </summary>
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

    /// <summary>
    /// Threshold for data. Will remove relative intensities below this relative to max intensity in each cluster
    /// </summary>
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

    /// <summary>
    /// Report multiple monoisotopic peaks
    /// </summary>
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
