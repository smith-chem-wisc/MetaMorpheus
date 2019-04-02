using System;
using System.Collections.ObjectModel;
using System.Linq;
using System.Windows.Controls;
using System.Windows.Threading;

namespace MetaMorpheusGUI
{
    class SearchModifications
    {
        public static DispatcherTimer Timer;
        public static bool FixedSearch;
        public static bool VariableSearch;
        public static bool GptmdSearch;
        
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
        public static void FilterTree(TextBox textbox, TreeView tree, ObservableCollection<ModTypeForTreeView> collection)
        {
            string key = textbox.Text.ToLower();
            if (string.IsNullOrEmpty(key))
            {
                tree.DataContext = collection; // shows full tree if nothing is searched
                return;
            }

            var modTypesWithMatchingMods = collection.Where(p => p.Children.Any(c => c.ModName.ToLower().Contains(key))); // parent of child mods that match key

            var modsThatMatchSearchString = new ObservableCollection<ModTypeForTreeView>(); // new collection containing expanded mod types that match key 

            foreach (ModTypeForTreeView modType in modTypesWithMatchingMods)
            {
                var textFilteredModType = new ModTypeForTreeView(modType.DisplayName, false);
                modsThatMatchSearchString.Add(textFilteredModType);
                textFilteredModType.Expanded = true;
                textFilteredModType.Use = modType.Use;

                var matchingChildren = modType.Children.Where(p => p.ModName.ToLower().Contains(key));
                foreach (ModForTreeView mod in matchingChildren)
                {
                    textFilteredModType.Children.Add(mod);
                }
            }

            tree.DataContext = modsThatMatchSearchString;
        }
    }
}