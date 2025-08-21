using Chemistry;
using Easy.Common.Extensions;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.Linq;
using MzLibUtil;
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
    public ObservableCollection<MassDifferenceAcceptorTypeViewModel> MassDiffAcceptorTypes { get; }
    public ObservableCollection<SelectableNotchViewModel> PredefinedNotches { get; } = new()
    {
        new SelectableNotchViewModel("Missed Mono", Chemistry.Constants.C13MinusC12)
    };

    public MassDifferenceAcceptorSelectionViewModel(MassDiffAcceptorType selectedType, string customText, double defaultCustomTolerance = 5)
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
        ToleranceValue = defaultCustomTolerance.ToString(CultureInfo.InvariantCulture);

        foreach (var notch in PredefinedNotches)
        {
            notch.PropertyChanged += (s, e) =>
            {
                if (e.PropertyName is nameof(SelectableNotchViewModel.IsSelected) or nameof(SelectableNotchViewModel.MaxPositiveFrequency) or nameof(SelectableNotchViewModel.MaxNegativeFrequency))
                        OnCustomFieldChanged();
            };
        }

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

    private MassDifferenceAcceptorTypeViewModel _selectedType;
    public MassDifferenceAcceptorTypeViewModel SelectedType
    {
        get => _selectedType;
        set
        {
            if (_selectedType != null) // only false during initialization when no cache is present. 
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
            {
                if (_cachedCustomMdac.IsNotNullOrEmpty())
                    CustomMdac = _cachedCustomMdac;
                OnCustomFieldChanged();
            }
        }
    }

    #region Custom 

    private bool _isManualCustomMdac = false;
    private string _cachedCustomMdac = string.Empty;
    private string _customMdac = null;
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

    public ObservableCollection<CustomMdacMode> CustomMdacModes { get; } = new(Enum.GetValues<CustomMdacMode>());
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
    private string _toleranceValue;
    public string ToleranceValue
    {
        get => _toleranceValue;
        set
        {
            _toleranceValue = value;
            OnPropertyChanged(nameof(ToleranceValue));
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
    private string _intervalRanges = "[-187,200]";
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

        string combinedDotMassShifts = DotMassShifts;
        if (CustomMode == CustomMdacMode.Notch)
        {
            // 1. Predefined notches
            var predefined = GetAllNotchCombinations()
                .Select(m => m.ToString("G6", System.Globalization.CultureInfo.InvariantCulture));

            // 2. User custom
            var user = Enumerable.Empty<string>();
            if (!string.IsNullOrWhiteSpace(DotMassShifts))
            {
                user = DotMassShifts
                    .Split(',', StringSplitOptions.RemoveEmptyEntries)
                    .Select(s => s.Trim())
                    .Where(s => double.TryParse(s, System.Globalization.NumberStyles.Any, System.Globalization.CultureInfo.InvariantCulture, out _));
            }

            // 3. Combine, deduplicate, sort
            var all = predefined.Concat(user)
                .Distinct()
                .OrderBy(x => double.Parse(x, System.Globalization.CultureInfo.InvariantCulture))
                .ToList();

            combinedDotMassShifts = string.Join(",", all);
        }

        return CustomMode switch
        {
            CustomMdacMode.Notch when ToleranceValue.IsNotNullOrEmpty() => $"{CustomName} dot {ToleranceValue} {SelectedToleranceType} {combinedDotMassShifts}",
            CustomMdacMode.Notch when ToleranceValue.IsNullOrEmpty() => $"{CustomName} dot {combinedDotMassShifts}",
            CustomMdacMode.Interval => $"{CustomName} interval {IntervalRanges}",
            CustomMdacMode.AroundZero => $"{CustomName} {SelectedToleranceType}AroundZero {ToleranceValue}",
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

                    const double tolerance = 1e-4;
                    if (split.Length >= 5)
                    {
                        ToleranceValue = split[2];
                        SelectedToleranceType = split[3].ToLowerInvariant();

                        // Parse the legacy mass shifts
                        var parsedShifts = split[4]
                            .Split(',', StringSplitOptions.RemoveEmptyEntries)
                            .Select(s => double.Parse(s, System.Globalization.CultureInfo.InvariantCulture))
                            .ToList();

                        // Try to match as many as possible using the notches
                        var notches = PredefinedNotches.ToList();
                        var maxPos = 5; 
                        var maxNeg = 5;

                        // Brute-force: try all combinations of frequencies for each notch (for a small number of notches and frequencies, this is feasible)
                        var bestMatch = new HashSet<double>();
                        int[] bestPos = new int[notches.Count];
                        int[] bestNeg = new int[notches.Count];

                        void TryCombinations(int notchIdx, int[] pos, int[] neg)
                        {
                            if (notchIdx == notches.Count)
                            {
                                // Generate all possible sums for this combination
                                var generated = new HashSet<double>();
                                var ranges = notches.Select((n, i) =>
                                    Enumerable.Range(-neg[i], pos[i] + neg[i] + 1).ToArray()).ToArray();

                                var combos = CartesianProduct(ranges);
                                foreach (var combo in combos)
                                {
                                    // Skip the all-zero combination (no mass shift)
                                    if (combo.All(x => x == 0))
                                        continue;

                                    double sum = 0;
                                    for (int i = 0; i < notches.Count; i++)
                                        sum += combo[i] * notches[i].MonoisotopicMass;
                                    generated.Add(sum);
                                }

                                var matched = parsedShifts.Where(s => generated.Any(g => Math.Abs(s - g) < tolerance)).ToHashSet();
                                if (matched.Count > bestMatch.Count)
                                {
                                    bestMatch = matched;
                                    Array.Copy(pos, bestPos, pos.Length);
                                    Array.Copy(neg, bestNeg, neg.Length);
                                }
                                return;
                            }
                            for (int p = 0; p <= maxPos; p++)
                            {
                                for (int n = 0; n <= maxNeg; n++)
                                {
                                    pos[notchIdx] = p;
                                    neg[notchIdx] = n;
                                    TryCombinations(notchIdx + 1, pos, neg);
                                }
                            }
                        }

                        TryCombinations(0, new int[notches.Count], new int[notches.Count]);

                        // Set the notches to the best found
                        for (int i = 0; i < notches.Count; i++)
                        {
                            notches[i].IsSelected = bestPos[i] > 0 || bestNeg[i] > 0;
                            notches[i].MaxPositiveFrequency = bestPos[i];
                            notches[i].MaxNegativeFrequency = bestNeg[i];
                        }

                        // Any unmatched shifts go into DotMassShifts
                        var unmatched = parsedShifts
                            .Where(s => !bestMatch.Any(m => Math.Abs(s - m) < tolerance))
                            .Select(s => s.ToString("G6", System.Globalization.CultureInfo.InvariantCulture))
                            .ToList();
                        DotMassShifts = string.Join(",", unmatched);
                    }
                    break;
                case "interval":
                    CustomMode = CustomMdacMode.Interval;
                    if (split.Length >= 3)
                        IntervalRanges = string.Join(' ', split.Skip(2));
                    break;
                case "ppmAroundZero":
                case "daltonsAroundZero":
                case "daAroundZero":
                    CustomMode = CustomMdacMode.AroundZero;
                    if (split.Length >= 3)
                        ToleranceValue = split[2];
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
        OnCustomFieldChanged();
    }

    private IEnumerable<double> GetAllNotchCombinations()
    {
        var selected = PredefinedNotches.Where(n => n.IsSelected).ToList();
        if (selected.Count == 0)
            yield return 0;
        else
        {
            // For each notch, create a list of possible counts (from -MaxNegative to +MaxPositive)
            var ranges = selected
                .Select(n => Enumerable.Range(-n.MaxNegativeFrequency, n.MaxPositiveFrequency + n.MaxNegativeFrequency + 1)
                    .ToArray())
                .ToArray();

            // Cartesian product of all ranges
            foreach (var combo in CartesianProduct(ranges))
            {
                double sum = 0;
                for (int i = 0; i < selected.Count; i++)
                    sum += combo[i] * selected[i].MonoisotopicMass;
                yield return sum;
            }
        }
    }

    // Helper for cartesian product
    private static IEnumerable<int[]> CartesianProduct(int[][] sequences)
    {
        int[] lengths = sequences.Select(s => s.Length).ToArray();
        int[] indices = new int[sequences.Length];
        while (true)
        {
            yield return indices.Select((v, i) => sequences[i][v]).ToArray();
            int k = sequences.Length - 1;
            while (k >= 0)
            {
                if (++indices[k] < lengths[k])
                    break;
                indices[k] = 0;
                k--;
            }
            if (k < 0)
                break;
        }
    }

    #endregion


    private MassDifferenceAcceptorTypeViewModel CreateModel(MassDiffAcceptorType type)
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
            MassDiffAcceptorType.PlusOrMinusThreeMM => "Basic search where the observed and theoretical precursor masses are allowed to disagree by +-1, +-2, or +-3 Da to allow for monoisotopic mass errors.",
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

        return new MassDifferenceAcceptorTypeViewModel
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
public class MassDifferenceAcceptorTypeViewModel : BaseViewModel, IEquatable<MassDiffAcceptorType>, IEquatable<MassDifferenceAcceptorTypeViewModel>
{
    private int _positiveMissedMonos;
    public int PositiveMissedMonos
    {
        get => _positiveMissedMonos;
        set
        {
            if (_positiveMissedMonos != value)
            {
                _positiveMissedMonos = value;
                OnPropertyChanged(nameof(PositiveMissedMonos));
            }
        }
    }

    private int _negativeMissedMonos;
    public int NegativeMissedMonos
    {
        get => _negativeMissedMonos;
        set
        {
            if (_negativeMissedMonos != value)
            {
                _negativeMissedMonos = value;
                OnPropertyChanged(nameof(NegativeMissedMonos));
            }
        }
    }

    public MassDiffAcceptorType Type { get; init; }
    public string Label { get; init; }
    public string ToolTip { get; init; }

    public bool Equals(MassDifferenceAcceptorTypeViewModel other)
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
        return Equals((MassDifferenceAcceptorTypeViewModel)obj);
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

public class SelectableNotchViewModel : BaseViewModel, IHasMass
{
    private bool _isSelected;
    private int _maxPositiveFrequency;
    private int _maxNegativeFrequency;
    public string Name { get; }
    public double MonoisotopicMass { get; }

    public int MaxPositiveFrequency
    {
        get => _maxPositiveFrequency;
        set { _maxPositiveFrequency = value; OnPropertyChanged(nameof(MaxPositiveFrequency)); }
    }

    public int MaxNegativeFrequency
    {
        get => _maxNegativeFrequency;
        set { _maxNegativeFrequency = value; OnPropertyChanged(nameof(MaxNegativeFrequency)); }
    }

    public bool IsSelected
    {
        get => _isSelected;
        set { _isSelected = value; OnPropertyChanged(nameof(IsSelected)); }
    }

    public SelectableNotchViewModel(string name, double mass, int maxFrequency = 1, int maxNegativeFrequency = 0)
    {
        Name = name;
        MonoisotopicMass = mass;
        MaxPositiveFrequency = maxFrequency;
        MaxNegativeFrequency = maxNegativeFrequency;
        IsSelected = false;
    }
}