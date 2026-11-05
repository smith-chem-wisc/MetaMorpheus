using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using Chemistry;
using EngineLayer;
using Omics.Fragmentation;
using Transcriptomics;
using Transcriptomics.Digestion;

namespace GuiFunctions;

/// <summary>
/// One choice in a terminal cap dropdown: a named cap, or the "custom formula" entry.
/// <see cref="Formula"/> is what is saved; it is null for the hydroxyl, which is what a null cap already means.
/// </summary>
/// <remarks>
/// Cap formulas are deltas to the sum of the residues, mzLib's convention: a residue is a nucleotide monophosphate
/// minus water, so it carries its own phosphate. That is why a 5' hydroxyl is "O-3P-1".
/// </remarks>
public class TerminalCapOption(string name, ChemicalFormula capFormula, string description, FragmentationTerminus end,
    bool isDefault = false, bool isCustom = false)
{
    public string Name { get; } = isDefault ? $"{name} (default)" : name;
    public ChemicalFormula CapFormula { get; } = capFormula;
    public string Formula => isDefault ? null : capFormula?.Formula;
    public bool IsDefault { get; } = isDefault;
    public bool IsCustom { get; } = isCustom;

    /// <summary>Monoisotopic shift relative to the hydroxyl at the same end.</summary>
    public double MassShift => capFormula == null ? 0 : capFormula.MonoisotopicMass -
        (end == FragmentationTerminus.FivePrime ? NucleicAcid.DefaultFivePrimeTerminus : NucleicAcid.DefaultThreePrimeTerminus).MonoisotopicMass;

    public string ToolTip => isCustom ? description
        : $"{description}\nFormula: {capFormula.Formula}   Shift vs. hydroxyl: {MassShift:+0.0000;-0.0000;0.0000} Da";

    public override string ToString() => Name;
}

/// <summary>
/// View model for choosing the 5' and 3' caps applied to RNA read from a database.
/// Only meaningful in RNA mode; a null cap means the mzLib default is used.
/// </summary>
public class RnaTerminalCapsViewModel : BaseViewModel
{
    public static readonly TerminalCapOption CustomOption = new("Custom formula...", null,
        "Type a chemical formula for the cap, e.g. H2O4P. Anything that cannot be read is ignored and the default is used.",
        FragmentationTerminus.None, isCustom: true);

    private TerminalCapOption _selectedFivePrimeOption;
    private TerminalCapOption _selectedThreePrimeOption;
    private string _fivePrimeCustomFormula;
    private string _threePrimeCustomFormula;

    public RnaTerminalCapsViewModel() : this(new CommonParameters())
    {
    }

    public RnaTerminalCapsViewModel(CommonParameters initialParams)
    {
        FivePrimeOptions = BuildOptions(FragmentationTerminus.FivePrime);
        ThreePrimeOptions = BuildOptions(FragmentationTerminus.ThreePrime);

        (_selectedFivePrimeOption, _fivePrimeCustomFormula) = Initialize(FivePrimeOptions, initialParams?.FivePrimeTerminus);
        (_selectedThreePrimeOption, _threePrimeCustomFormula) = Initialize(ThreePrimeOptions, initialParams?.ThreePrimeTerminus);
    }

    // Shifts checked against published values (5'/3' phosphate +79.97, 3' cyclic phosphate +61.96) and HPO3 arithmetic.
    private static ObservableCollection<TerminalCapOption> BuildOptions(FragmentationTerminus end)
    {
        TerminalCapOption Cap(string name, string formula, string description, bool isDefault = false)
            => new(name, ChemicalFormula.ParseFormula(formula), description, end, isDefault);

        var options = end == FragmentationTerminus.FivePrime
            ?
            [
                Cap("Hydroxyl", NucleicAcid.DefaultFivePrimeTerminus.Formula, "5'-OH. The default: the first residue's phosphate is removed.", true),
                Cap("Monophosphate", "H", "5'-phosphate, as left by many nucleases and by kinase treatment."),
                Cap("Diphosphate", "H2O3P", "5'-diphosphate (5'-pp)."),
                Cap("Triphosphate", "H3O6P2", "5'-triphosphate (5'-ppp), as on in vitro transcripts without capping.")
            ]
            : new List<TerminalCapOption>
            {
                Cap("Hydroxyl", NucleicAcid.DefaultThreePrimeTerminus.Formula, "3'-OH. The default for a full-length oligonucleotide.", true),
                Cap("Phosphate", Rnase.DefaultThreePrimeTerminus.ThisChemicalFormula.Formula, "3'-phosphate on the full-length RNA's 3' end only. Digestion products already end in a 3'-phosphate."),
                Cap("Cyclic phosphate", "O3P", "2',3'-cyclic phosphate on the full-length RNA's 3' end only. Interior digestion products are not changed.")
            };
        options.Add(CustomOption);
        return new ObservableCollection<TerminalCapOption>(options);
    }

    private static (TerminalCapOption, string) Initialize(IList<TerminalCapOption> options, ChemicalFormula formula)
    {
        if (formula == null)
            return (options[0], string.Empty); // the hydroxyl is always first

        // match on the formula itself: its text can be re-ordered ("H1O1" comes back as "HO")
        var match = options.FirstOrDefault(o => !o.IsCustom && formula.Equals(o.CapFormula));
        return match != null ? (match, string.Empty) : (CustomOption, formula.Formula);
    }

    public ObservableCollection<TerminalCapOption> FivePrimeOptions { get; }
    public ObservableCollection<TerminalCapOption> ThreePrimeOptions { get; }

    public TerminalCapOption SelectedFivePrimeOption
    {
        get => _selectedFivePrimeOption;
        set
        {
            _selectedFivePrimeOption = value;
            OnPropertyChanged(nameof(SelectedFivePrimeOption));
            OnPropertyChanged(nameof(IsFivePrimeCustom));
        }
    }

    public TerminalCapOption SelectedThreePrimeOption
    {
        get => _selectedThreePrimeOption;
        set
        {
            _selectedThreePrimeOption = value;
            OnPropertyChanged(nameof(SelectedThreePrimeOption));
            OnPropertyChanged(nameof(IsThreePrimeCustom));
        }
    }

    public string FivePrimeCustomFormula
    {
        get => _fivePrimeCustomFormula;
        set
        {
            _fivePrimeCustomFormula = value;
            OnPropertyChanged(nameof(FivePrimeCustomFormula));
        }
    }

    public string ThreePrimeCustomFormula
    {
        get => _threePrimeCustomFormula;
        set
        {
            _threePrimeCustomFormula = value;
            OnPropertyChanged(nameof(ThreePrimeCustomFormula));
        }
    }

    public bool IsFivePrimeCustom => SelectedFivePrimeOption?.IsCustom == true;
    public bool IsThreePrimeCustom => SelectedThreePrimeOption?.IsCustom == true;

    /// <summary>The chosen 5' cap formula, or null for the mzLib default.</summary>
    public string FivePrimeTerminus => Resolve(SelectedFivePrimeOption, FivePrimeCustomFormula);

    /// <summary>The chosen 3' cap formula, or null for the mzLib default.</summary>
    public string ThreePrimeTerminus => Resolve(SelectedThreePrimeOption, ThreePrimeCustomFormula);

    private static string Resolve(TerminalCapOption option, string customFormula)
    {
        if (option == null)
            return null;
        if (!option.IsCustom)
            return option.Formula;
        return string.IsNullOrWhiteSpace(customFormula) ? null : customFormula.Trim();
    }

    /// <summary>
    /// Checks that any custom formula can be parsed. Returns false and a message for the first invalid cap.
    /// </summary>
    public bool Validate(out string error)
    {
        error = null;
        foreach (var (label, formula) in new[] { ("5'", FivePrimeTerminus), ("3'", ThreePrimeTerminus) })
        {
            if (formula == null)
                continue;
            try
            {
                ChemicalFormula.ParseFormula(formula);
            }
            catch (Exception)
            {
                error = $"The {label} cap '{formula}' is not a valid chemical formula (for example: {(label == "5'" ? "H or H2O3P" : "H2O4P or O3P")}). Leave it blank to use the default.";
                return false;
            }
        }
        return true;
    }
}
