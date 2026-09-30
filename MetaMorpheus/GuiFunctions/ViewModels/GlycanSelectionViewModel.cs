using EngineLayer;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.Linq;
using System.Text;

namespace GuiFunctions;

/// <summary>
/// The glycan tree in the Glyco Search task window: one group per glycan database, every glycan in
/// it a checkable row, so a search can be narrowed to individual glycans instead of a whole database.
/// </summary>
/// <remarks>
/// The glycans are handed in rather than read from <c>GlobalVariables</c> so a test can supply its
/// own, matching how the other task-window view models take the task's parameter objects.
///
/// Reuses <see cref="ModTypeForTreeViewModel"/> and <see cref="ModForTreeViewModel"/> unchanged: a
/// Glycan is a Modification, so the existing nodes already fit and no glycan-specific node type is
/// needed. What this adds is the wiring those nodes lack -- see <see cref="OnGlycanUseChanged"/>.
/// </remarks>
public class GlycanSelectionViewModel : BaseViewModel
{
    private bool _syncingCheckState;
    private string _summary;
    private HashSet<string> _activeDatabases; // null until the window says, which means every database
    private string _filterKey = string.Empty;

    /// <summary>
    /// Master group -> the copy of it <see cref="ModTreeFilter.Filter"/> made for the tree while the
    /// search box is filtering. The copies share their rows with the master groups, but not their
    /// check state, so both have to be settled when a row changes.
    /// </summary>
    private readonly Dictionary<ModTypeForTreeViewModel, ModTypeForTreeViewModel> _filteredCopies =
        new Dictionary<ModTypeForTreeViewModel, ModTypeForTreeViewModel>();

    /// <summary>
    /// One node per glycan database, in the order given: the master collection the read-back reads,
    /// never a filtered view of it. Holds every database, including those the search will not use.
    /// </summary>
    public ObservableCollection<ModTypeForTreeViewModel> Databases { get; } = new ObservableCollection<ModTypeForTreeViewModel>();

    /// <summary>
    /// What the tree shows, and so its ItemsSource: the groups for the databases the search will use
    /// (see <see cref="SetActiveDatabases"/>), narrowed by <see cref="Filter"/> when the search box
    /// holds text. The same instance throughout, so the tree is bound once.
    /// </summary>
    public ObservableCollection<ModTypeForTreeViewModel> Displayed { get; } = new ObservableCollection<ModTypeForTreeViewModel>();

    /// <summary>
    /// The master groups for the databases the search will use.
    /// </summary>
    private IEnumerable<ModTypeForTreeViewModel> ActiveDatabases =>
        _activeDatabases == null ? Databases : Databases.Where(db => _activeDatabases.Contains(db.DisplayName));

    /// <summary>
    /// What is checked, in words, for the line beside the search box.
    /// </summary>
    public string Summary
    {
        get => _summary;
        private set
        {
            _summary = value;
            OnPropertyChanged(nameof(Summary));
        }
    }

    /// <param name="glycansByDatabase">
    /// Database file name -> the glycans it holds. <c>GlobalVariables.OGlycansByDatabase</c>
    /// concatenated with <c>NGlycansByDatabase</c> is what the task window passes.
    /// </param>
    /// <param name="selectedGlycans">
    /// A previously saved selection as (database file name, glycan IdWithMotif) pairs. Null or empty
    /// means nothing is checked, which in turn means the whole database will be searched.
    /// </param>
    public GlycanSelectionViewModel(
        IEnumerable<KeyValuePair<string, List<Glycan>>> glycansByDatabase,
        IEnumerable<(string, string)> selectedGlycans = null)
    {
        foreach (var database in glycansByDatabase ?? Enumerable.Empty<KeyValuePair<string, List<Glycan>>>())
        {
            // DisplayName is the bare file name on purpose: it is half of the persisted key that
            // ToSelectedGlycans projects and the engine matches on, so decorating it -- with a count,
            // say -- would break both.
            var theDatabase = new ModTypeForTreeViewModel(database.Key, false);
            Databases.Add(theDatabase);

            foreach (var glycan in database.Value ?? new List<Glycan>())
            {
                AddGlycan(theDatabase, GlycanToolTip(glycan), glycan.IdWithMotif, GlycanLabel(glycan), bad: false);
            }

            SettleCheckState(theDatabase);
        }

        Restore(selectedGlycans);
        RefreshDisplayed();
        UpdateSummary();
    }

    /// <summary>
    /// Limits the tree, the summary and the saved selection to the databases the search will use: the
    /// one each glycan-database combo box names, for the search type chosen.
    /// </summary>
    /// <remarks>
    /// The engine narrows only the databases it searches, so a tick in any other group would be saved,
    /// counted in the summary, and then ignored. Ticks in a hidden group are kept, not cleared, so
    /// choosing that database again brings them back.
    /// </remarks>
    public void SetActiveDatabases(IEnumerable<string> databaseFileNames)
    {
        _activeDatabases = new HashSet<string>((databaseFileNames ?? Enumerable.Empty<string>()).Where(name => name != null));
        RefreshDisplayed();
        UpdateSummary();
    }

    /// <summary>
    /// Narrows <see cref="Displayed"/> to the rows whose label contains <paramref name="key"/>,
    /// case-insensitively. Empty shows every row again.
    /// </summary>
    /// <remarks>
    /// Matches the row label rather than the identifier, so "HexNAc", "2222" and "H5N4A2" all work.
    /// Owned here rather than by the window so the filtered copies can be kept in step with their rows
    /// -- see <see cref="OnGlycanUseChanged"/>.
    /// </remarks>
    public void Filter(string key)
    {
        _filterKey = key ?? string.Empty;
        RefreshDisplayed();
    }

    private void RefreshDisplayed()
    {
        _filteredCopies.Clear();
        Displayed.Clear();

        if (string.IsNullOrEmpty(_filterKey))
        {
            foreach (var database in ActiveDatabases)
            {
                Displayed.Add(database);
            }
            return;
        }

        foreach (var copy in ModTreeFilter.Filter(ActiveDatabases, _filterKey, row => row.DisplayName))
        {
            // Filter only returns a group holding at least one match, and its rows are the master's own.
            _filteredCopies[copy.Children[0].Parent] = copy;
            Displayed.Add(copy);
        }

        // Filter copies the master group's check state, but the copy's checkbox acts only on the rows it
        // shows, so it has to read from those. Guarded because settling a group cascades down to its rows.
        _syncingCheckState = true;
        try
        {
            foreach (var copy in _filteredCopies.Values)
            {
                SettleCheckState(copy);
            }
        }
        finally
        {
            _syncingCheckState = false;
        }
    }

    /// <summary>
    /// Re-checks the glycans a saved task named, and marks in red any it names that are no longer
    /// there -- the same way the modification trees surface a mod this install does not have.
    /// </summary>
    private void Restore(IEnumerable<(string, string)> selectedGlycans)
    {
        foreach (var (databaseName, glycanId) in selectedGlycans ?? Enumerable.Empty<(string, string)>())
        {
            var theDatabase = Databases.FirstOrDefault(b => b.DisplayName.Equals(databaseName));
            if (theDatabase == null)
            {
                // The database itself is gone, so the group has to be invented to hang the row on.
                theDatabase = new ModTypeForTreeViewModel(databaseName, true);
                Databases.Add(theDatabase);
                AddGlycan(theDatabase, UnknownGlycan, glycanId, glycanId, bad: true);
            }
            else
            {
                var theGlycan = theDatabase.Children.FirstOrDefault(b => b.ModName.Equals(glycanId));
                if (theGlycan != null)
                {
                    theGlycan.Use = true;
                }
                else
                {
                    AddGlycan(theDatabase, UnknownGlycan, glycanId, glycanId, bad: true);
                }
            }

            SettleCheckState(theDatabase);
        }
    }

    private void AddGlycan(ModTypeForTreeViewModel database, string toolTip, string modName, string displayName, bool bad)
    {
        var node = new ModForTreeViewModel(toolTip, bad, modName, bad, database, displayName);
        node.PropertyChanged += OnGlycanUseChanged;
        database.Children.Add(node);
    }

    /// <summary>
    /// Keeps the database checkbox and <see cref="Summary"/> in step with the glycans underneath them.
    /// </summary>
    /// <remarks>
    /// <see cref="ModForTreeViewModel.Use"/> notifies only itself; nothing tells the parent a child
    /// changed, and <c>VerifyCheckState</c> is a manual pull the owner has to make. The modification
    /// trees only ever make it while restoring a saved task, which is why a group there still reads
    /// as fully checked right after you uncheck something inside it.
    ///
    /// The re-entrancy guard is required, not defensive: VerifyCheckState assigns the parent's Use,
    /// and a non-null assignment cascades back down to every child, each of which raises
    /// PropertyChanged again. Checking the last box in a group would otherwise recurse forever.
    /// </remarks>
    private void OnGlycanUseChanged(object sender, PropertyChangedEventArgs e)
    {
        if (_syncingCheckState || e.PropertyName != nameof(ModForTreeViewModel.Use))
        {
            return;
        }

        _syncingCheckState = true;
        try
        {
            if (sender is ModForTreeViewModel glycan)
            {
                SettleCheckState(glycan.Parent);

                // While filtering, the checkbox on screen is the copy, not the master group.
                if (glycan.Parent != null && _filteredCopies.TryGetValue(glycan.Parent, out var copy))
                {
                    SettleCheckState(copy);
                }
            }

            UpdateSummary();
        }
        finally
        {
            _syncingCheckState = false;
        }
    }

    /// <summary>
    /// Brings one database checkbox into agreement with its children.
    /// </summary>
    /// <remarks>
    /// The empty case has to be said outright. <see cref="ModTypeForTreeViewModel"/> takes <c>bad</c>
    /// as its second constructor argument, not <c>use</c>, so Use is never initialized and a bool?
    /// defaults to null -- which is the indeterminate box. VerifyCheckState over no children leaves
    /// it there, so a database holding no glycans would be the only one that looked partly selected.
    /// </remarks>
    private static void SettleCheckState(ModTypeForTreeViewModel database)
    {
        if (database == null)
        {
            return;
        }

        if (database.Children.Count == 0)
        {
            database.Use = false;
        }
        else
        {
            database.VerifyCheckState();
        }
    }

    private void UpdateSummary()
    {
        var active = ActiveDatabases.ToList();
        int selected = active.Sum(db => db.Children.Count(c => c.Use));
        int total = active.Sum(db => db.Children.Count);
        int databases = active.Count(db => db.Children.Any(c => c.Use));

        // "entries", not "glycans": a glycan is listed once per attachment site, so the 12 lines of
        // OGlycan.gdb are 24 rows here. Calling them glycans invites the reader to compare the number
        // with the file and find it wrong.
        Summary = selected == 0
            ? string.Format(CultureInfo.InvariantCulture,
                "   none checked " + EmDash + " the whole selected database will be searched ({0} entries available)", total)
            : string.Format(CultureInfo.InvariantCulture,
                "   {0} of {1} entries checked, across {2} database{3}", selected, total, databases, databases == 1 ? "" : "s");
    }

    /// <summary>
    /// The checked glycans as (database file name, glycan IdWithMotif) pairs, for
    /// <c>GlycoSearchParameters.SelectedGlycans</c>. No pairs for a database means the whole database.
    /// </summary>
    /// <remarks>
    /// Only leaves are read, and always from <see cref="Databases"/> rather than from whatever the
    /// tree is currently displaying -- so a selection made while the search box is filtered is still
    /// returned in full.
    ///
    /// Two kinds of group write nothing. A fully checked group is saved as the whole database, not as
    /// the entries it holds today: a list would silently leave out anything added to the file later,
    /// and the custom databases exist to be added to. A group for a database the search will not use
    /// is left out, because the engine would ignore it.
    ///
    /// Keyed by composition string rather than by Glycan.GlyId, which is a positional index into the
    /// loaded array: a saved index would quietly point at a different glycan if the .gdb were edited.
    /// </remarks>
    public List<(string, string)> ToSelectedGlycans()
    {
        return ActiveDatabases
            .Where(database => database.Use != true)
            .SelectMany(database => database.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.ModName)))
            .ToList();
    }

    /// <summary>
    /// The row label: identifier, expanded composition, and mass.
    /// </summary>
    /// <remarks>
    /// IdWithMotif alone is a composition code ("H5N4A2 on N"), unreadable unless you already know
    /// the letters and unsearchable in the terms people use -- typing "HexNAc" would match nothing.
    /// Spelling the monosaccharides out and appending the mass makes both work, because the tree
    /// searches this label.
    /// </remarks>
    public static string GlycanLabel(Glycan glycan)
    {
        var expanded = ExpandComposition(glycan.Composition);
        var mass = glycan.Mass / 1E5; // Glycan.Mass is the monoisotopic mass scaled by 1e5

        return string.IsNullOrEmpty(expanded)
            ? string.Format(CultureInfo.InvariantCulture, "{0}   {1}   {2:F2} Da", glycan.IdWithMotif, EmDash, mass)
            : string.Format(CultureInfo.InvariantCulture, "{0}   {1}   {2}   {1}   {3:F2} Da", glycan.IdWithMotif, EmDash, expanded, mass);
    }

    /// <summary>
    /// Written as an escape rather than the character itself: the literal was once saved through a
    /// Windows-1252 round trip and came out as the three characters U+00E2 U+20AC U+201D.
    /// </summary>
    private const string EmDash = "\u2014";

    public static string GlycanToolTip(Glycan glycan)
    {
        var lines = new List<string>
        {
            "Composition: " + glycan.Composition,
            "Expanded:    " + ExpandComposition(glycan.Composition),
            string.Format(CultureInfo.InvariantCulture, "Mass:        {0:F5} Da", glycan.Mass / 1E5),
            "Type:        " + glycan.Type,
        };

        // Only structure-format databases carry a structure string; composition databases leave it null.
        if (!string.IsNullOrEmpty(glycan.Struc))
        {
            lines.Add("Structure:   " + glycan.Struc);
        }

        return string.Join(Environment.NewLine, lines);
    }

    /// <summary>
    /// Turns a composition code such as "H5N4A2" into "Hex(5)HexNAc(4)NeuAc(2)".
    /// </summary>
    /// <remarks>
    /// Built from Glycan.NameCharDic so custom monosaccharides registered at start-up are spelled out
    /// too, rather than silently falling back to their letter. A code with no count means one, and an
    /// unrecognized code is passed through as itself.
    /// </remarks>
    public static string ExpandComposition(string composition)
    {
        if (string.IsNullOrWhiteSpace(composition))
        {
            return string.Empty;
        }

        // NameCharDic holds aliases as well as canonical names -- "dHex" and "Fuc" share a code and a
        // slot -- so an inversion has to choose, and must not choose by dictionary order. Ordering by
        // slot then ordinally by name is deterministic and picks the canonical name, which is also the
        // spelling the shipped databases use on disk.
        var nameByCode = Glycan.NameCharDic
            .GroupBy(entry => entry.Value.Item1)
            .ToDictionary(
                group => group.Key,
                group => group.OrderBy(entry => entry.Value.Item2)
                              .ThenBy(entry => entry.Key, StringComparer.Ordinal)
                              .First().Key);

        var expanded = new StringBuilder();
        int i = 0;
        while (i < composition.Length)
        {
            char code = composition[i++];
            int start = i;
            while (i < composition.Length && char.IsDigit(composition[i]))
            {
                i++;
            }

            var count = i > start ? composition.Substring(start, i - start) : "1";
            expanded.Append(nameByCode.TryGetValue(code, out var name) ? name : code.ToString())
                    .Append('(').Append(count).Append(')');
        }

        return expanded.ToString();
    }

    private const string UnknownGlycan = "UNKNOWN GLYCAN!";
}

[ExcludeFromCodeCoverage] // For design time gui display only
public class GlycanSelectionModel : GlycanSelectionViewModel
{
    public static GlycanSelectionModel Instance => new GlycanSelectionModel();

    public GlycanSelectionModel() : base(new Dictionary<string, List<Glycan>>())
    {
    }
}
