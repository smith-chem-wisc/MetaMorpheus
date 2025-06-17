using System;
using System.Collections.ObjectModel;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using TaskLayer;

namespace GuiFunctions;
public enum CustomMdacMode
{
    Notch,
    Interval, 
    AroundZero, // <-- single mode for both ppm and da
}

/// <summary>
/// Controls the selection of different mass diff acceptors in the GUI and Custom Mdac handling. 
/// </summary>
public class MassDifferenceAcceptorSelectionViewModel : BaseViewModel
{
    public ObservableCollection<MassDifferenceAcceptorTypeModel> MassDiffAcceptorTypes { get; }
    public MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType selectedType, string customText) : base()
    {
        // Almost every piece of this constructor is order dependent. Be careful when changing. 
        var types = Enum.GetValues<MassDiffAcceptorType>().ToList();
        if (types.Count > 1)
        {
            var last = types[^1];
            types.RemoveAt(types.Count - 1);
            types.Insert(4, last);
        }
        var models = types.Select(CreateModel);
        MassDiffAcceptorTypes = new(models);
        SelectedType = MassDiffAcceptorTypes.FirstOrDefault(m => m.Type == selectedType) ?? MassDiffAcceptorTypes.First();

        // Default tolerance types
        ToleranceTypes = new[] { "ppm", "da" };
        SelectedToleranceType = ToleranceTypes[0];

        // Parse legacy string if present
        if (!string.IsNullOrWhiteSpace(customText))
        {
            ParseLegacyCustomMdac(customText);
        }
        else
        {
            CustomMdac = BuildCustomMdac();
        }
    }

    private MassDifferenceAcceptorTypeModel _selectedType;
    public MassDifferenceAcceptorTypeModel SelectedType
    {
        get => _selectedType;
        set
        {
            if (_selectedType != null) // only true during initialization when no cache is present. 
            {
                switch (_selectedType.Type)
                {
                    // if moving away from a string input type, cache it
                    case MassDiffAcceptorType.Custom:
                        _cachedCustomMdac = CustomMdac;
                        CustomMdac = string.Empty;
                        break;
                }
            }

            _selectedType = value;
            OnPropertyChanged(nameof(SelectedType));

            // if moving to string input type, populate mdac with cached string
            if (_selectedType.Type == MassDiffAcceptorType.Custom)
                CustomMdac = _cachedCustomMdac;
        }
    }

    #region Custom 

    private bool _isManualCustomMdac = false;
    private string _cachedCustomMdac = string.Empty;
    private string _customMdac = string.Empty;
    public string CustomMdac
    {
        get => _isManualCustomMdac ? _customMdac : BuildCustomMdac();
        set
        {
            if (_customMdac != value)
            {
                _customMdac = value;
                _isManualCustomMdac = true;
                ParseLegacyCustomMdac(value);
                OnPropertyChanged(nameof(CustomMdac));
            }
        }
    }

    public string[] ToleranceTypes { get; } 
    private string _selectedToleranceType;
    public string SelectedToleranceType
    {
        get => _selectedToleranceType;
        set
        {
            _selectedToleranceType = value;
            OnPropertyChanged(nameof(SelectedToleranceType));
            OnCustomFieldChanged();
        }
    }

    public ObservableCollection<CustomMdacMode> CustomMdacMods { get; set; } = new(Enum.GetValues<CustomMdacMode>());
    private CustomMdacMode _customMode;
    public CustomMdacMode CustomMode
    {
        get => _customMode;
        set
        {
            _customMode = value;
            OnPropertyChanged(nameof(CustomMode));
            OnCustomFieldChanged();
        }
    }

    // Common name for all custom types
    private string _customName = "Custom";
    public string CustomName
    {
        get => _customName;
        set
        {
            _customName = value;
            OnPropertyChanged(nameof(CustomName));
            OnCustomFieldChanged();
        }
    }

    // Notch mode
    private string _dotTolerance = "5";
    public string DotTolerance
    {
        get => _dotTolerance;
        set
        {
            _dotTolerance = value;
            OnPropertyChanged(nameof(DotTolerance));
            OnCustomFieldChanged();
        }
    }

    private string _dotMassShifts = "0";
    public string DotMassShifts
    {
        get => _dotMassShifts;
        set
        {
            _dotMassShifts = value;
            OnPropertyChanged(nameof(DotMassShifts));
            OnCustomFieldChanged();
        }
    }

    // Interval mode
    private string _intervalRanges;
    public string IntervalRanges
    {
        get => _intervalRanges;
        set
        {
            _intervalRanges = value;
            OnPropertyChanged(nameof(IntervalRanges));
            OnCustomFieldChanged();
        }
    }

    // AroundZero mode - ppm or dalton tolerance value. 
    private string _aroundZeroValue;
    public string AroundZeroValue
    {
        get => _aroundZeroValue;
        set
        {
            _aroundZeroValue = value;
            OnPropertyChanged(nameof(AroundZeroValue));
            OnCustomFieldChanged();
        }
    }

    // Helper to update CustomMdac when any field changes
    private void OnCustomFieldChanged()
    {
        _isManualCustomMdac = false;
        OnPropertyChanged(nameof(CustomMdac));
    }


    private string BuildCustomMdac()
    {
        if (SelectedType?.Type != MassDiffAcceptorType.Custom)
            return string.Empty;

        return CustomMode switch
        {
            CustomMdacMode.Notch => $"{CustomName} dot {DotTolerance} {SelectedToleranceType} {DotMassShifts}",
            CustomMdacMode.Interval => $"{CustomName} interval {IntervalRanges}",
            CustomMdacMode.AroundZero => $"{CustomName} {SelectedToleranceType}AroundZero {AroundZeroValue}",
            _ => _customMdac
        };
    }

    /// <summary>
    /// Parses a legacy customMdac string and populates the fields for the new GUI.
    /// </summary>
    private void ParseLegacyCustomMdac(string text)
    {
        try
        {
            var split = text.Split(' ', StringSplitOptions.RemoveEmptyEntries);
            if (split.Length < 2)
                return;

            CustomName = split[0];

            switch (split[1])
            {
                case "dot":
                    CustomMode = CustomMdacMode.Notch;
                    if (split.Length >= 5)
                    {
                        DotTolerance = split[2];
                        SelectedToleranceType = split[3].ToLowerInvariant();
                        DotMassShifts = split[4];
                    }
                    break;
                case "interval":
                    CustomMode = CustomMdacMode.Interval;
                    if (split.Length >= 3)
                        IntervalRanges = string.Join(' ', split.Skip(2));
                    break;
                case "ppmAroundZero":
                case "daltonsAroundZero":
                    CustomMode = CustomMdacMode.AroundZero;
                    if (split.Length >= 3)
                        AroundZeroValue = split[2];
                    SelectedToleranceType = split[1].ToLower().StartsWith("ppm") ? "ppm" : "da";
                    break;
                case "OpenSearch":
                    // Switch to the Open type and do not set CustomMode
                    SelectedType = MassDiffAcceptorTypes.First(p => p.Type == MassDiffAcceptorType.Open);
                    break;
                default:
                    // For unknown/unsupported modes, just store the string for advanced/legacy use
                    _customMdac = text;
                    break;
            }
        }
        catch
        {
            // fallback: just store the string
            _customMdac = text;
        }
    }

    #endregion



    private MassDifferenceAcceptorTypeModel CreateModel(MassDiffAcceptorType type)
    {
        string label = type switch
        {
            MassDiffAcceptorType.Exact => "Exact",
            MassDiffAcceptorType.OneMM => "1 Missed Monoisotopic Peak",
            MassDiffAcceptorType.TwoMM => "1 or 2 Missed Monoisotopic Peaks",
            MassDiffAcceptorType.ThreeMM => "1, 2, or 3 Missed Monoisotopic Peaks",
            MassDiffAcceptorType.PlusOrMinusThreeMM => "+- 3 Missed Monoisotopic Peaks",
            MassDiffAcceptorType.ModOpen => "-187 and Up",
            MassDiffAcceptorType.Open => "Accept all",
            MassDiffAcceptorType.Custom => "Custom",
            _ => throw new NotImplementedException($"No model implemented for type: {type}"),
        };

        string toolTip = type switch
        {
            MassDiffAcceptorType.Exact => "Basic search where the observed and theoretical precursor masses must be equal (~0 Da precursor mass-difference). This search type assumes that there are no monoisotopic errors.",
            MassDiffAcceptorType.OneMM => "Basic search where the observed and theoretical precursor masses are allowed to disagree by 1 Da to allow for a 1 Da monoisotopic mass error.",
            MassDiffAcceptorType.TwoMM => "Basic search where the observed and theoretical precursor masses are allowed to disagree by 1 or 2 Da to allow for a 1 or 2 Da monoisotopic mass error.",
            MassDiffAcceptorType.ThreeMM => "Basic search where the observed and theoretical precursor masses are allowed to disagree by 1, 2, or 3 Da to allow for a 1, 2, or 3 Da monoisotopic mass error.",
            MassDiffAcceptorType.PlusOrMinusThreeMM => "Basic search where the observed and theoretical precursor masses are allowed to disagree by +-1, +-2, or +-3 Da in to allow for monoisotopic mass errors.",
            MassDiffAcceptorType.ModOpen => "An \"open-mass\" search that allows mass-differences between observed and theoretical precursor masses of -187 Da to infinity (observed can be infinitely more massive than the theoretical).\r\nThe purpose of this search type is to detect mass-differences corresponding to PTMs, amino acid variants, sample handling artifacts, etc.\r\nPlease use \"Modern Search\" mode when using this search type.",
            MassDiffAcceptorType.Open => "An \"open-mass\" search that allows mass-differences between observed and theoretical precursor masses of -infinity to infinity. \r\nThe purpose of this search type is to detect mass-differences corresponding to PTMs, amino acid variants, sample handling artifacts, etc. \r\nPlease use \"Modern Search\" mode when using this search type.",
            MassDiffAcceptorType.Custom => "A custom mass difference acceptor may be specified in multiple ways: \r\n* To accept a custom (other than the interval corresponding to the precursor tolerance) interval around zero daltons, specify a custom name, followed by \"ppmAroundZero\" or \"daltonsAroundZero\", followed by the numeric value corresponding to the interval width. \r\nExamples: \r\n\t* CustomPpmInterval ppmAroundZero 5 \r\n\t* CustomDaltonInterval daltonsAroundZero 2.1 \r\n* To accept a variety of pre-specified mass differences, use a custom name, followed by \"dot\", followed by a custom bin width, followed by comma separated acceptable mass differences. \r\nExamples: \r\n\t* CustomMissedIsotopePeaks dot 5 ppm 0,1.0029,2.0052 \r\n\t* CustomOxidationAllowed dot 0.1 da 0,16 \r\n* To accept mass differences in pre-specified dalton intervals, use a custom name, followed by \"interval\", followed by comma separated mass intervals in brackets. \r\nExample: \r\n\t* CustomPositiveIntervalAcceptror interval [0,200]",
            _ => throw new NotImplementedException($"No model implemented for type: {type}"),
        };

        int positiveMissedMonos = type switch
        {
            MassDiffAcceptorType.Exact => 0,
            MassDiffAcceptorType.OneMM => 1,
            MassDiffAcceptorType.TwoMM => 2,
            MassDiffAcceptorType.ThreeMM => 3,
            MassDiffAcceptorType.PlusOrMinusThreeMM => 3,
            MassDiffAcceptorType.ModOpen => 0,
            MassDiffAcceptorType.Open => 0,
            MassDiffAcceptorType.Custom => 0,
            _ => throw new NotImplementedException($"No model implemented for type: {type}"),
        };

        int negativeMissedMonos = type switch
        {
            MassDiffAcceptorType.Exact => 0,
            MassDiffAcceptorType.OneMM => 0,
            MassDiffAcceptorType.TwoMM => 0,
            MassDiffAcceptorType.ThreeMM => 0,
            MassDiffAcceptorType.PlusOrMinusThreeMM => 3,
            MassDiffAcceptorType.ModOpen => 0,
            MassDiffAcceptorType.Open => 0,
            MassDiffAcceptorType.Custom => 0,
            _ => throw new NotImplementedException($"No model implemented for type: {type}"),
        };

        return new MassDifferenceAcceptorTypeModel
        {
            Type = type,
            Label = label,
            ToolTip = toolTip,
            PositiveMissedMonos = positiveMissedMonos,
            NegativeMissedMonos = negativeMissedMonos
        };
    }
}

/// <summary>
/// A single MassDiff acceptor type and all relevant GUI display information
/// </summary>
public class MassDifferenceAcceptorTypeModel : IEquatable<MassDiffAcceptorType>, IEquatable<MassDifferenceAcceptorTypeModel>
{
    public int PositiveMissedMonos { get; init; }
    public int NegativeMissedMonos { get; init; }
    public MassDiffAcceptorType Type { get; init; }
    public string Label { get; init; }
    public string ToolTip { get; init; }

    public bool Equals(MassDifferenceAcceptorTypeModel other)
    {
        return Type == other!.Type;
    }

    public bool Equals(MassDiffAcceptorType other)
    {
        return Type == other;
    }

    public override bool Equals(object obj)
    {
        if (obj is null) return false;
        if (ReferenceEquals(this, obj)) return true;
        if (obj.GetType() != GetType()) return false;
        return Equals((MassDifferenceAcceptorTypeModel)obj);
    }

    public override int GetHashCode()
    {
        return (int)Type;
    }
}

[ExcludeFromCodeCoverage] // For design time gui display only
public class MassDifferenceAcceptorSelectionModel : MassDifferenceAcceptorSelectionViewModel
{
    public static MassDifferenceAcceptorSelectionModel Instance => new MassDifferenceAcceptorSelectionModel();
    public MassDifferenceAcceptorSelectionModel() : base(MassDiffAcceptorType.TwoMM, "")
    {

    }
}