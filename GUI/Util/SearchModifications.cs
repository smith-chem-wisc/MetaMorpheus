using System;
using System.Collections.ObjectModel;
using System.Linq;
using System.Windows.Controls;
using System.Windows.Threading;

namespace MetaMorpheusGUI
{
    class SearchModifications
    {
        public DispatcherTimer Timer;
        public bool FixedSearch { get; set; } // true if fixed mods search box used
        public bool VarSearch { get; set; } // true if var mods search box used

        public SearchModifications()
        {
            Timer = new DispatcherTimer();
            Timer.Interval = TimeSpan.FromMilliseconds(300);
        }

        // starts timer to keep track of user keystrokes
        public void SetTimer()
        {
            Timer.Stop(); // Resets the timer
            Timer.Start();
        }

        // filters and expands tree according to user mod search
        public static void FilterTree(TextBox textbox, TreeView tree, ObservableCollection<ModTypeForTreeView> collection)
        {
            string key = textbox.Text.ToLower();
            if (key == "")
            {
                tree.DataContext = collection; // shows full tree if nothing is searched
                return;
            }

            var parentModsMatch = collection.Where(p => p.DisplayName.ToLower().Contains(key)); // parent mod types that match key
            var childModsMatch = collection.Where(p => p.Children.Any(c => c.DisplayName.ToLower().Contains(key))); // parent of child mods that match key

            if (parentModsMatch.Count() == 0 && childModsMatch.Count() != 0)
            {
                var allModsMatch = new ObservableCollection<ModTypeForTreeView>(); // new collection containing expanded mod types that match key                  
                foreach (ModTypeForTreeView parent in childModsMatch.ToList())
                {
                    var newParent = new ModTypeForTreeView(parent.DisplayName, false);
                    allModsMatch.Add(newParent);
                    newParent.Use = false;
                    newParent.Expanded = true;

                    var children = parent.Children.Where(y => y.DisplayName.ToLower().Contains(key));
                    children.ToList().ForEach(child => newParent.Children.Add(child));
                }
                parentModsMatch = allModsMatch;
            }
            tree.DataContext = parentModsMatch;
        }
    }
}