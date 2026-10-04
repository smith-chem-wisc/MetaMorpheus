using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Windows.Input;
using Easy.Common.Extensions;
using EngineLayer;
using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using TaskLayer;
using Transcriptomics;
using Transcriptomics.Digestion;
using GuiFunctions.Util;

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
    private bool _modsCanSuppressBaseLossFragments;

    public FragmentationParamsViewModel(CommonParameters initialParams, SearchParameters searchParams)
    {
        // Initialize from provided parameters
        _generateComplementaryIons = initialParams.AddCompIons;
        if (initialParams.FragmentationParameters is not null)
        {
            _generateMIon = initialParams.FragmentationParameters.GenerateMIon;
            _modsCanSuppressBaseLossFragments = initialParams.FragmentationParameters is RnaFragmentationParams
            {
                ModificationsCanSuppressBaseLossIons: true
            };
        }
        else
        {
            _modsCanSuppressBaseLossFragments = false;
            _generateMIon = GuiGlobalParamsViewModel.Instance.IsRnaMode;
        }

        _maxFragmentMassDa = searchParams.MaxFragmentSize;
        _minInternalIonLength = searchParams.MinAllowedInternalFragmentLength;
        _useInternalIons = searchParams.MinAllowedInternalFragmentLength > 0;

        // Initialize fragment ion termini
        var terminus = initialParams.DigestionParams.FragmentationTerminus;
        _leftSideFragmentIons = terminus == FragmentationTerminus.FivePrime || terminus == FragmentationTerminus.N || terminus == FragmentationTerminus.Both;
        _rightSideFragmentIons = terminus == FragmentationTerminus.ThreePrime || terminus == FragmentationTerminus.C || terminus == FragmentationTerminus.Both;

        // Load available M-ion losses (built-in + custom)
        var availableMIonLosses = LoadAvailableMIonLosses(initialParams);
        AvailableMIonLosses = new ObservableCollection<MIonLossViewModel>(availableMIonLosses);

        AddAllLossesCommand = new RelayCommand(AddAllLosses);
        RemoveAllLossesCommand = new RelayCommand(RemoveAllLosses);
    }

    private List<MIonLossViewModel> LoadAvailableMIonLosses(CommonParameters initialParams)
    {
        var allLosses = new List<MIonLossViewModel>();
        
        // Determine which losses were previously selected
        var selectedLosses = initialParams.FragmentationParameters?.MIonLosses ?? new List<MIonLoss>();

        if (GuiGlobalParamsViewModel.Instance.IsRnaMode)
        {
            // Ensure static constructor has run to populate MIonLoss.AllMIonLosses
            _ = new RnaFragmentationParams();
            
            // Add built-in RNA losses
            foreach (var loss in MIonLoss.AllMIonLosses.Values)
            {
                bool isSelected = selectedLosses.Any(l => l.Annotation == loss.Annotation);
                allLosses.Add(new MIonLossViewModel(loss, isSelected, false));
            }
        }
        else
        {
            // Add built-in peptide losses (if any exist)
            if (MIonLoss.AllMIonLosses.TryGetValue("-H2O", out var waterLoss))
            {
                bool isSelected = selectedLosses.Any(l => l.Annotation == waterLoss.Annotation);
                allLosses.Add(new MIonLossViewModel(waterLoss, isSelected, true));
            }
            if (MIonLoss.AllMIonLosses.TryGetValue("-NH3", out var ammoniaLoss))
            {
                bool isSelected = selectedLosses.Any(l => l.Annotation == ammoniaLoss.Annotation);
                allLosses.Add(new MIonLossViewModel(ammoniaLoss, isSelected, true));
            }
        }

        // Add custom losses
        var customLosses = CustomMIonLossManager.LoadCustomMIonLosses()
            .Where(l => l.IsApplicableToCurrentMode());

        foreach (var customLoss in customLosses)
        {
            var mIonLoss = customLoss.ToMIonLoss();
            bool isSelected = selectedLosses.Any(l => l.Annotation == customLoss.Annotation);
            allLosses.Add(new MIonLossViewModel(mIonLoss, isSelected, true));
        }

        return allLosses;
    }

    public ICommand AddAllLossesCommand { get; set; }
    public ICommand RemoveAllLossesCommand { get; set; }

    public bool GenerateMIon
    {
        get => _generateMIon;
        set
        {
            _generateMIon = value;
            OnPropertyChanged(nameof(GenerateMIon));
        }
    }

    public bool ModsCanSuppressBaseLossFragments
    {
        get => _modsCanSuppressBaseLossFragments;
        set
        {
            _modsCanSuppressBaseLossFragments = value;
            OnPropertyChanged(nameof(ModsCanSuppressBaseLossFragments));
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

    public ObservableCollection<MIonLossViewModel> AvailableMIonLosses { get; protected set; }

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
                ModificationsCanSuppressBaseLossIons = ModsCanSuppressBaseLossFragments,
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

    /// <summary>
    /// Reloads the M-Ion losses from the cache/file
    /// </summary>
    public void ReloadMIonLosses()
    {
        // Force reload by clearing cache
        CustomMIonLossManager.ClearCache();
        
        // Preserve current terminus state to avoid resetting user preferences
        var terminus = LeftSideFragmentIons && RightSideFragmentIons ? FragmentationTerminus.Both
            : LeftSideFragmentIons ? FragmentationTerminus.N
            : RightSideFragmentIons ? FragmentationTerminus.C
            : FragmentationTerminus.Both;

        // Build CommonParameters with correct digeston type to avoid CommonParameters
        // overriding the passed fragmentationParams with a protein default
        IDigestionParams digParams = GuiGlobalParamsViewModel.Instance?.IsRnaMode == true
            ? new RnaDigestionParams(fragmentationTerminus: terminus)
            : new DigestionParams(fragmentationTerminus: terminus);

        var commonParams = new CommonParameters(
            digestionParams: digParams,
            fragmentationParams: ToFragmentationParams()
        );

        // LoadAvailableMIonLosses already restores selections from FragmentationParameters.MIonLosses
        // No need to manually restore selections here
        AvailableMIonLosses.Clear();
        foreach (var loss in LoadAvailableMIonLosses(commonParams))
            AvailableMIonLosses.Add(loss);
        
        OnPropertyChanged(nameof(AvailableMIonLosses));
    }

    private void AddAllLosses() => AvailableMIonLosses.ForEach(p => p.IsSelected = true);
    private void RemoveAllLosses() => AvailableMIonLosses.ForEach(p => p.IsSelected = false);
}

/// <summary>
/// View model wrapper for MIonLoss to enable selection in the GUI
/// </summary>
public class MIonLossViewModel : BaseViewModel
{
    public MIonLossViewModel(MIonLoss mIonLoss, bool isSelected = false, bool isCustom = false)
    {
        MIonLoss = mIonLoss;
        IsSelected = isSelected;
        IsCustom = isCustom;
    }

    public MIonLoss MIonLoss { get; }

    public string Name => MIonLoss.Name;
    public string Annotation => MIonLoss.Annotation;
    public double MassShift => MIonLoss.MonoisotopicMass;
    public string ChemicalFormula => MIonLoss.ThisChemicalFormula.Formula;
    public bool IsCustom { get; }

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
        
        // Add some sample losses for design-time viewing
        AvailableMIonLosses = new ObservableCollection<MIonLossViewModel>
        {
            new MIonLossViewModel(new MIonLoss("Water Loss", "-H2O", Chemistry.ChemicalFormula.ParseFormula("H2O1")), true, true),
            new MIonLossViewModel(new MIonLoss("Ammonia Loss", "-NH3", Chemistry.ChemicalFormula.ParseFormula("H3N1")), false, true),
            new MIonLossViewModel(new MIonLoss("Phosphate Loss", "-P", Chemistry.ChemicalFormula.ParseFormula("H1P1O3")), true, false),
        };
    }
}
