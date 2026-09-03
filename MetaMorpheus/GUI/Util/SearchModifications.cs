using System;
using System.Collections.ObjectModel;
using System.Linq;
using System.Windows.Controls;
using System.Windows.Threading;
using GuiFunctions;

namespace MetaMorpheusGUI
{
    class SearchModifications
    {
        public static DispatcherTimer Timer;
        public static bool FixedSearch;
        public static bool VariableSearch;
        public static bool GptmdSearch;
        public static bool GlycanSearch;
        
        public static void SetUpModSearchBoxes()
        {
            Timer = new DispatcherTimer();
            Timer.Interval = TimeSpan.FromMilliseconds(300);
        }

        // starts timer to keep track of user keystrokes
        public static void SetTimer()
        {
            // Reset the timer
            Timer.Stop();
            Timer.Start();
        }

        // filters and expands tree according to user mod search
        public static void FilterTree(TextBox textbox, TreeView tree, ObservableCollection<ModTypeForTreeViewModel> collection)
        {
            FilterTree(textbox, tree, collection, m => m.ModName);
        }

        /// <summary>
        /// As <see cref="FilterTree(TextBox,TreeView,ObservableCollection{ModTypeForTreeViewModel})"/>, but
        /// searching text of the caller's choosing instead of only <see cref="ModForTreeViewModel.ModName"/>.
        /// </summary>
        /// <remarks>
        /// Modification names are English words, so matching the identifier is enough -- "ribosyl" finds
        /// "ADP-ribosylation on S". A glycan's identifier is a composition code, so matching it alone means
        /// "HexNAc" finds nothing even in a database that spells HexNAc(1) on disk. This lets the glycan tree
        /// search the expanded name and mass as well, without changing what the other search boxes do.
        /// </remarks>
        public static void FilterTree(TextBox textbox, TreeView tree, ObservableCollection<ModTypeForTreeViewModel> collection,
            Func<ModForTreeViewModel, string> searchText)
        {
            string key = textbox.Text.ToLower();
            if (string.IsNullOrEmpty(key))
            {
                tree.DataContext = collection; // shows full tree if nothing is searched
                return;
            }

            var modTypesWithMatchingMods = collection.Where(p => p.Children.Any(c => Matches(c, searchText, key))); // parent of child mods that match key

            var modsThatMatchSearchString = new ObservableCollection<ModTypeForTreeViewModel>(); // new collection containing expanded mod types that match key 

            foreach (ModTypeForTreeViewModel modType in modTypesWithMatchingMods)
            {
                var textFilteredModType = new ModTypeForTreeViewModel(modType.DisplayName, false);
                modsThatMatchSearchString.Add(textFilteredModType);
                textFilteredModType.Expanded = true;
                textFilteredModType.Use = modType.Use;

                var matchingChildren = modType.Children.Where(p => Matches(p, searchText, key));
                foreach (ModForTreeViewModel mod in matchingChildren)
                {
                    textFilteredModType.Children.Add(mod);
                }
            }

            tree.DataContext = modsThatMatchSearchString;
        }

        private static bool Matches(ModForTreeViewModel mod, Func<ModForTreeViewModel, string> searchText, string key)
        {
            var text = searchText(mod);
            return text != null && text.ToLower().Contains(key);
        }
    }
}