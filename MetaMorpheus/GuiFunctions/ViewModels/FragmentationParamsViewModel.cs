using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using EngineLayer;
using Omics.Fragmentation;
using TaskLayer;
using Transcriptomics;

namespace GuiFunctions;

/// <summary>
/// View model for managing fragmentation parameters in the GUI.
/// Contains properties specific to RNA fragmentation when in RNA mode.
/// </summary>
public class FragmentationParamsViewModel : BaseViewModel
{
    private bool _generateMIon;
    private bool _generateComplementaryIons;
    private bool _leftSideFragmentIons;
    private bool _rightSideFragmentIons;
    private double _maxFragmentMassDa;
    private int _minInternalIonLength;
    private bool _useInternalIons;

    public FragmentationParamsViewModel(CommonParameters initialParams, SearchParameters searchParams)
    {
        // Initialize from provided parameters
        _generateMIon = initialParams.FragmentationParameters?.GenerateMIon ?? GuiGlobalParamsViewModel.Instance.IsRnaMode;
        _generateComplementaryIons = initialParams.AddCompIons;
        _maxFragmentMassDa = searchParams.MaxFragmentSize;
        _minInternalIonLength = searchParams.MinAllowedInternalFragmentLength;
        _useInternalIons = searchParams.MinAllowedInternalFragmentLength > 0;

        // Initialize fragment ion termini
        var terminus = initialParams.DigestionParams.FragmentationTerminus;
        _leftSideFragmentIons = terminus == FragmentationTerminus.FivePrime || terminus == FragmentationTerminus.N || terminus == FragmentationTerminus.Both;
        _rightSideFragmentIons = terminus == FragmentationTerminus.ThreePrime || terminus == FragmentationTerminus.C || terminus == FragmentationTerminus.Both;

        // Initialize available M-ion losses from the static collection
        if (GuiGlobalParamsViewModel.Instance.IsRnaMode)
        {
            _ = new RnaFragmentationParams(); // Ensure static constructor has run to populate MIonLoss.AllMIonLosses
            AvailableMIonLosses = new ObservableCollection<MIonLossViewModel>(
                MIonLoss.AllMIonLosses.Values.Select(m => new MIonLossViewModel(m, (initialParams.FragmentationParameters as RnaFragmentationParams)!.MIonLosses.Contains(m))));
        }
        else 
        {
            AvailableMIonLosses = new ObservableCollection<MIonLossViewModel>();
        }
    }

    public bool GenerateMIon
    {
        get => _generateMIon;
        set
        {
            _generateMIon = value;
            OnPropertyChanged(nameof(GenerateMIon));
        }
    }

    public bool GenerateComplementaryIons
    {
        get => _generateComplementaryIons;
        set
        {
            _generateComplementaryIons = value;
            OnPropertyChanged(nameof(GenerateComplementaryIons));
        }
    }

    public bool LeftSideFragmentIons
    {
        get => _leftSideFragmentIons;
        set
        {
            _leftSideFragmentIons = value;
            OnPropertyChanged(nameof(LeftSideFragmentIons));
        }
    }

    public bool RightSideFragmentIons
    {
        get => _rightSideFragmentIons;
        set
        {
            _rightSideFragmentIons = value;
            OnPropertyChanged(nameof(RightSideFragmentIons));
        }
    }

    public double MaxFragmentMassDa
    {
        get => _maxFragmentMassDa;
        set
        {
            _maxFragmentMassDa = value;
            OnPropertyChanged(nameof(MaxFragmentMassDa));
        }
    }

    public int MinInternalIonLength
    {
        get => _minInternalIonLength;
        set
        {
            _minInternalIonLength = value;
            OnPropertyChanged(nameof(MinInternalIonLength));
        }
    }

    public bool GenerateInternalIons
    {
        get => _useInternalIons;
        set
        {
            _useInternalIons = value;
            if (value && MinInternalIonLength <= 0)
            {
                // Set a default minimum length for internal ions if enabling them without a valid length
                MinInternalIonLength = 4; 
            }
            if (!value && MinInternalIonLength > 0)
            {
                // Reset minimum length if disabling internal ions
                MinInternalIonLength = 0;
            }
            OnPropertyChanged(nameof(GenerateInternalIons));
        }
    }

    public ObservableCollection<MIonLossViewModel> AvailableMIonLosses { get; private set; }

    /// <summary>
    /// Gets the selected M-ion losses as a list
    /// </summary>
    public List<MIonLoss> GetSelectedMIonLosses()
    {
        return AvailableMIonLosses.Where(m => m.IsSelected).Select(m => m.MIonLoss).ToList();
    }

    /// <summary>
    /// Converts the view model to FragmentationParams object
    /// </summary>
    public IFragmentationParams ToFragmentationParams()
    {
        if (GuiGlobalParamsViewModel.Instance.IsRnaMode)
        {
            return new RnaFragmentationParams
            {
                GenerateMIon = GenerateMIon,
                MIonLosses = GetSelectedMIonLosses()
            };
        }
        else
        {
            return new FragmentationParams
            {
                GenerateMIon = GenerateMIon,
                MIonLosses = GetSelectedMIonLosses()
            };
        }
    }
}

/// <summary>
/// View model wrapper for MIonLoss to enable selection in the GUI
/// </summary>
public class MIonLossViewModel : BaseViewModel
{
    public MIonLossViewModel(MIonLoss mIonLoss, bool isSelected = false)
    {
        MIonLoss = mIonLoss;
        IsSelected = isSelected;
    }

    public MIonLoss MIonLoss { get; }

    public string Name => MIonLoss.Name;
    public string Annotation => MIonLoss.Annotation;

    private bool _isSelected;
    public bool IsSelected
    {
        get => _isSelected;
        set
        {
            _isSelected = value;
            OnPropertyChanged(nameof(IsSelected));
        }
    }
}

[ExcludeFromCodeCoverage] // Design-time model for XAML designer
public class FragmentationParamsDesignModel : FragmentationParamsViewModel
{
    public static FragmentationParamsDesignModel Instance => new();

    public FragmentationParamsDesignModel() : base(new(), new())
    {
        GenerateMIon = true;
    }
}
