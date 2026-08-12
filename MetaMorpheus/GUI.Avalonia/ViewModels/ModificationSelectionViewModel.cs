using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using CommunityToolkit.Mvvm.ComponentModel;
using EngineLayer;
using Omics.Modifications;

namespace MetaMorpheus.Avalonia.ViewModels;

/// <summary>
/// Chooses fixed and variable modifications.
///
/// A modification's identity in CommonParameters is the pair (ModificationType, IdWithMotif) - not the
/// Modification object - because that is what round-trips through TOML. So selection is tracked as
/// those pairs, and matched back against GlobalVariables.AllModsKnown for display.
///
/// GlobalVariables loads over 4,000 modifications, which is far too many for a flat list, so they are
/// grouped by type and filterable. That grouping is also how the WPF window presents them.
/// </summary>
internal sealed partial class ModificationSelectionViewModel : ObservableObject
{
    private readonly List<ModificationChoice> _all;

    [ObservableProperty] private string _filter = string.Empty;

    public ObservableCollection<ModificationGroup> Groups { get; } = new();

    public ModificationSelectionViewModel(
        IEnumerable<(string, string)> fixedMods, IEnumerable<(string, string)> variableMods)
    {
        var chosenFixed = new HashSet<(string, string)>(fixedMods ?? Enumerable.Empty<(string, string)>());
        var chosenVariable = new HashSet<(string, string)>(variableMods ?? Enumerable.Empty<(string, string)>());

        _all = GlobalVariables.AllModsKnown
            .Where(m => m.ValidModification && m.ModificationType is not null && m.IdWithMotif is not null)
            .GroupBy(m => (m.ModificationType, m.IdWithMotif))
            .Select(g => g.First())
            .Select(m => new ModificationChoice(m)
            {
                IsFixed = chosenFixed.Contains((m.ModificationType, m.IdWithMotif)),
                IsVariable = chosenVariable.Contains((m.ModificationType, m.IdWithMotif)),
            })
            .OrderBy(c => c.ModificationType, StringComparer.Ordinal)
            .ThenBy(c => c.Name, StringComparer.Ordinal)
            .ToList();

        Rebuild();
    }

    public int TotalCount => _all.Count;

    public int SelectedFixedCount => _all.Count(c => c.IsFixed);

    public int SelectedVariableCount => _all.Count(c => c.IsVariable);

    partial void OnFilterChanged(string value) => Rebuild();

    /// <summary>
    /// Rebuilds the grouped view for the current filter. Anything already ticked stays visible even if
    /// it does not match, so a filter cannot hide a selection the user has made and then forget it.
    /// </summary>
    private void Rebuild()
    {
        string needle = Filter?.Trim() ?? string.Empty;

        IEnumerable<ModificationChoice> visible = _all.Where(c =>
            c.IsFixed || c.IsVariable
            || needle.Length == 0
            || c.Name.Contains(needle, StringComparison.OrdinalIgnoreCase)
            || c.ModificationType.Contains(needle, StringComparison.OrdinalIgnoreCase));

        Groups.Clear();
        foreach (var group in visible.GroupBy(c => c.ModificationType).OrderBy(g => g.Key, StringComparer.Ordinal))
        {
            Groups.Add(new ModificationGroup(group.Key, group.ToList()));
        }
    }

    public IReadOnlyList<(string, string)> FixedSelection =>
        _all.Where(c => c.IsFixed).Select(c => (c.ModificationType, c.Name)).ToList();

    public IReadOnlyList<(string, string)> VariableSelection =>
        _all.Where(c => c.IsVariable).Select(c => (c.ModificationType, c.Name)).ToList();

    /// <summary>
    /// A modification cannot sensibly be both fixed and variable, so the two are mutually exclusive -
    /// the WPF window allows it, which produces a search that treats the residue inconsistently.
    /// </summary>
    public void SetFixed(ModificationChoice choice, bool isFixed)
    {
        choice.IsFixed = isFixed;
        if (isFixed)
        {
            choice.IsVariable = false;
        }
        OnPropertyChanged(nameof(SelectedFixedCount));
        OnPropertyChanged(nameof(SelectedVariableCount));
    }

    public void SetVariable(ModificationChoice choice, bool isVariable)
    {
        choice.IsVariable = isVariable;
        if (isVariable)
        {
            choice.IsFixed = false;
        }
        OnPropertyChanged(nameof(SelectedFixedCount));
        OnPropertyChanged(nameof(SelectedVariableCount));
    }

    /// <summary>Restores the defaults MetaMorpheus uses when no modifications are specified.</summary>
    public void ResetToDefaults()
    {
        foreach (ModificationChoice choice in _all)
        {
            choice.IsFixed = choice.ModificationType == "Common Fixed"
                && (choice.Name == "Carbamidomethyl on C" || choice.Name == "Carbamidomethyl on U");
            choice.IsVariable = choice.ModificationType == "Common Variable" && choice.Name == "Oxidation on M";
        }
        Rebuild();
        OnPropertyChanged(nameof(SelectedFixedCount));
        OnPropertyChanged(nameof(SelectedVariableCount));
    }
}

/// <summary>One selectable modification. Name is IdWithMotif, which is the half stored in TOML.</summary>
internal sealed partial class ModificationChoice : ObservableObject
{
    [ObservableProperty] private bool _isFixed;
    [ObservableProperty] private bool _isVariable;

    public ModificationChoice(Modification modification)
    {
        ModificationType = modification.ModificationType;
        Name = modification.IdWithMotif;
        MonoisotopicMass = modification.MonoisotopicMass;
    }

    public string ModificationType { get; }
    public string Name { get; }
    public double? MonoisotopicMass { get; }

    public string MassDisplay => MonoisotopicMass.HasValue ? $"{MonoisotopicMass.Value:F4}" : string.Empty;
}

internal sealed class ModificationGroup
{
    public ModificationGroup(string modificationType, IReadOnlyList<ModificationChoice> choices)
    {
        ModificationType = modificationType;
        Choices = choices;
    }

    public string ModificationType { get; }
    public IReadOnlyList<ModificationChoice> Choices { get; }
    public string Header => $"{ModificationType} ({Choices.Count})";
}
