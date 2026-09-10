using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;

namespace GuiFunctions;

/// <summary>
/// Narrows a modification/glycan tree to the nodes matching what the user typed.
/// </summary>
/// <remarks>
/// Lives here rather than beside the search boxes in the GUI project so it can be tested: the
/// caller keeps the WPF part (read the TextBox, assign the TreeView's DataContext) and this keeps
/// the part with the behaviour in it.
/// </remarks>
public static class ModTreeFilter
{
    /// <summary>
    /// The text a row is matched on when the caller does not say. Modification identifiers are
    /// English words, so matching the identifier is enough for the modification trees.
    /// </summary>
    public static string DefaultSearchText(ModForTreeViewModel mod) => mod.ModName;

    /// <param name="collection">The master collection. It is not modified.</param>
    /// <param name="key">What the user typed. Matching is case-insensitive substring.</param>
    /// <param name="searchText">
    /// The text of a row to match against, defaulting to <see cref="DefaultSearchText"/>. The glycan
    /// tree passes the display label instead, because a glycan identifier is a composition code --
    /// matching only that means "HexNAc" finds nothing, even in a database that spells HexNAc(1) on
    /// disk but parses to N1.
    /// </param>
    /// <returns>
    /// Groups that contain at least one match, holding only their matching children and expanded.
    /// </returns>
    public static ObservableCollection<ModTypeForTreeViewModel> Filter(
        IEnumerable<ModTypeForTreeViewModel> collection,
        string key,
        Func<ModForTreeViewModel, string> searchText = null)
    {
        searchText ??= DefaultSearchText;
        var lowered = (key ?? string.Empty).ToLower();

        var matching = new ObservableCollection<ModTypeForTreeViewModel>();

        foreach (var modType in collection.Where(p => p.Children.Any(c => Matches(c, searchText, lowered))))
        {
            var filtered = new ModTypeForTreeViewModel(modType.DisplayName, false);
            matching.Add(filtered);

            // Order matters. Use is copied while Children is still empty, because the Use setter
            // cascades to every child -- writing it after the children were added would overwrite
            // the real check state, and the children below are the SAME instances as the master
            // collection's, which is what makes filtering non-destructive.
            filtered.Expanded = true;
            filtered.Use = modType.Use;

            foreach (var mod in modType.Children.Where(p => Matches(p, searchText, lowered)))
            {
                filtered.Children.Add(mod);
            }
        }

        return matching;
    }

    private static bool Matches(ModForTreeViewModel mod, Func<ModForTreeViewModel, string> searchText, string loweredKey)
    {
        var text = searchText(mod);
        return text != null && text.ToLower().Contains(loweredKey);
    }
}
