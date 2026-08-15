using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using CommunityToolkit.Mvvm.ComponentModel;
using EngineLayer;
using Omics.Digestion;
using Omics.Modifications;

namespace MetaMorpheusAvalonia.ViewModels;

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

    /// <summary>
    /// <paramref name="isRna"/> selects which modification set is offered. GlobalVariables keeps two -
    /// AllModsKnown and AllRnaModsKnown - and an RNA task's modifications are only in the second, so
    /// offering AllModsKnown for one means none of its modifications can be chosen at all.
    /// </summary>
    public ModificationSelectionViewModel(
        IEnumerable<(string, string)> fixedMods,
        IEnumerable<(string, string)> variableMods,
        IEnumerable<(string, string)> gptmdMods = null,
        bool isRna = false)
    {
        IsRna = isRna;
        var chosenFixed = new HashSet<(string, string)>(fixedMods ?? Enumerable.Empty<(string, string)>());
        var chosenVariable = new HashSet<(string, string)>(variableMods ?? Enumerable.Empty<(string, string)>());
        var chosenGptmd = new HashSet<(string, string)>(gptmdMods ?? Enumerable.Empty<(string, string)>());

        // One entry per (ModificationType, IdWithMotif), because that pair is the identity TOML stores -
        // two modifications sharing it cannot be selected independently at all. On the shipped
        // databases this collapses 84 groups, all in the glycan sets, and CollapsedDuplicates
        // confirms every one of them is a byte-identical duplicate rather than a distinct variant.
        var byIdentity = (isRna ? GlobalVariables.AllRnaModsKnown : GlobalVariables.AllModsKnown)
            .Where(m => m.ValidModification && m.ModificationType is not null && m.IdWithMotif is not null)
            .GroupBy(m => (m.ModificationType, m.IdWithMotif))
            .ToList();

        CollapsedDuplicates = byIdentity.Count(g => g.Count() > 1);
        DistinctModificationsCollapsed = byIdentity.Count(HoldsDistinctModifications);

        _all = byIdentity
            .Select(g => g.First())
            .Select(m => new ModificationChoice(this, m)
            {
                IsFixed = chosenFixed.Contains((m.ModificationType, m.IdWithMotif)),
                IsVariable = chosenVariable.Contains((m.ModificationType, m.IdWithMotif)),
                IsGptmd = chosenGptmd.Contains((m.ModificationType, m.IdWithMotif)),
            })
            .OrderBy(c => c.ModificationType, StringComparer.Ordinal)
            .ThenBy(c => c.Name, StringComparer.Ordinal)
            .ToList();

        Rebuild();
    }

    /// <summary>
    /// True when entries sharing a TOML identity are not interchangeable - which would mean collapsing
    /// them hides a modification the user cannot otherwise reach. Zero on the shipped databases; a
    /// test asserts that, so a database change that broke the assumption would be caught rather than
    /// silently narrowing the choices.
    /// </summary>
    private static bool HoldsDistinctModifications(IGrouping<(string, string), Modification> group) =>
        group.Select(m => m.LocationRestriction).Distinct().Count() > 1
        || group.Select(m => m.MonoisotopicMass).Distinct().Count() > 1
        || group.Select(m => m.Target?.ToString()).Distinct().Count() > 1;

    /// <summary>Which of GlobalVariables' two modification sets is being offered.</summary>
    public bool IsRna { get; }

    /// <summary>How many TOML identities matched more than one entry in the offered set.</summary>
    public int CollapsedDuplicates { get; }

    /// <summary>How many of those held entries that genuinely differ. Expected to be zero.</summary>
    public int DistinctModificationsCollapsed { get; }

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
    /// GPTMD's ListOfModsGptmd - the one setting users actually change on a GPTMD task. Independent of
    /// fixed and variable, since GPTMD is asking which modifications to go looking for.
    /// </summary>
    public IReadOnlyList<(string, string)> GptmdSelection =>
        _all.Where(c => c.IsGptmd).Select(c => (c.ModificationType, c.Name)).ToList();

    public int SelectedGptmdCount => _all.Count(c => c.IsGptmd);

    /// <summary>
    /// A modification cannot sensibly be both fixed and variable, so the two are mutually exclusive.
    /// The WPF window allows it - SearchTaskWindow.xaml.cs harvests two independent trees with no
    /// cross-check - which produces a search that treats the residue inconsistently. This is a
    /// deliberate divergence from the WPF GUI, called out in the PR description.
    ///
    /// The rule lives in ModificationChoice's property setters rather than here, so that a XAML author
    /// binding IsChecked straight to IsFixed cannot get around it.
    /// </summary>
    public void SetFixed(ModificationChoice choice, bool isFixed) => choice.IsFixed = isFixed;

    public void SetVariable(ModificationChoice choice, bool isVariable) => choice.IsVariable = isVariable;

    /// <summary>Raised by ModificationChoice whenever a selection changes, however it was set.</summary>
    internal void NoteSelectionChanged()
    {
        OnPropertyChanged(nameof(SelectedFixedCount));
        OnPropertyChanged(nameof(SelectedVariableCount));
        OnPropertyChanged(nameof(SelectedGptmdCount));
    }

    /// <summary>
    /// Restores the defaults the engine itself would apply. Read off a CommonParameters rather than
    /// written out here, so they cannot drift from CommonParameters.cs - which also picks different
    /// defaults for RNA, hence the digestion parameters being passed through.
    /// </summary>
    public void ResetToDefaults(IDigestionParams digestionParams = null)
    {
        var defaults = new CommonParameters(digestionParams: digestionParams);
        var defaultFixed = new HashSet<(string, string)>(defaults.ListOfModsFixed);
        var defaultVariable = new HashSet<(string, string)>(defaults.ListOfModsVariable);

        foreach (ModificationChoice choice in _all)
        {
            var identity = (choice.ModificationType, choice.Name);
            choice.SetWithoutExclusion(
                isFixed: defaultFixed.Contains(identity),
                isVariable: defaultVariable.Contains(identity));
        }

        Rebuild();
        NoteSelectionChanged();
    }
}

/// <summary>One selectable modification. Name is IdWithMotif, which is the half stored in TOML.</summary>
internal sealed partial class ModificationChoice : ObservableObject
{
    private readonly ModificationSelectionViewModel _owner;
    private bool _applyingDefaults;

    [ObservableProperty] private bool _isFixed;
    [ObservableProperty] private bool _isVariable;
    [ObservableProperty] private bool _isGptmd;

    public ModificationChoice(ModificationSelectionViewModel owner, Modification modification)
    {
        _owner = owner;
        ModificationType = modification.ModificationType;
        Name = modification.IdWithMotif;
        MonoisotopicMass = modification.MonoisotopicMass;
    }

    // Exclusivity is enforced here rather than in the parent so that it holds however the property is
    // set - including <CheckBox IsChecked="{Binding IsFixed}"/>, the binding a XAML author reaches for
    // first, which bypassed the parent's methods entirely.
    partial void OnIsFixedChanged(bool value)
    {
        if (value && !_applyingDefaults)
        {
            IsVariable = false;
        }
        _owner?.NoteSelectionChanged();
    }

    partial void OnIsVariableChanged(bool value)
    {
        if (value && !_applyingDefaults)
        {
            IsFixed = false;
        }
        _owner?.NoteSelectionChanged();
    }

    // GPTMD is not exclusive with either: it asks which modifications to search for, not how to apply
    // one that is already known to be there.
    partial void OnIsGptmdChanged(bool value) => _owner?.NoteSelectionChanged();

    /// <summary>Sets both flags as a unit, for ResetToDefaults - which is already consistent.</summary>
    internal void SetWithoutExclusion(bool isFixed, bool isVariable)
    {
        _applyingDefaults = true;
        try
        {
            IsFixed = isFixed;
            IsVariable = isVariable;
        }
        finally
        {
            _applyingDefaults = false;
        }
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
