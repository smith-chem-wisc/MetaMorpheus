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
            }
            else
            {
                var parentMod = collection.Where(p => p.DisplayName.ToLower().Contains(key));
                if (parentMod.Count() == 0) // no parent mod type matches key
                {
                    var modMatches = new ObservableCollection<ModTypeForTreeView>(); // new collection containing mod types that match key
                    var parentChildMatch = collection.Where(p => p.Children.Any(c => c.DisplayName.ToLower().Contains(key))); // parent mod types of children that matches key

                    if (parentChildMatch.Count() != 0)
                    {
                        foreach (ModTypeForTreeView parent in parentChildMatch.ToList())
                        {
                            var newParent = new ModTypeForTreeView(parent.DisplayName, false);
                            modMatches.Add(newParent);
                            newParent.Use = false;
                            newParent.Expanded = true;
                            
                            var children = parent.Children.Where(y => y.DisplayName.ToLower().Contains(key));
                            foreach (ModForTreeView child in children)
                            {
                                newParent.Children.Add(child);
                            }
                        }
                        parentMod = modMatches;
                    } 
                }
                tree.DataContext = parentMod;
            }
        }
    }
}