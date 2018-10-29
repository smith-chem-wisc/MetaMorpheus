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
        public bool TimerCreated { get; set; }
        public bool FixedSearch { get; set; }
        public bool VarSearch { get; set; }

        public SearchModifications()
        {
            Timer = new DispatcherTimer();
            Timer.Interval = TimeSpan.FromMilliseconds(300);
            TimerCreated = false;
        }

        public void SetTimer()
        {
            Timer.Stop(); // Resets the timer
            Timer.Start();
        }
        
        public void FilterTree(TextBox textbox, TreeView tree, ObservableCollection<ModTypeForTreeView> collection)
        {
            if (textbox.Text == "")
            {
                tree.DataContext = collection; // shows full tree if nothing is searched
            }
            else
            {
                var parent = collection.Where(d => d.DisplayName.ToLower().Contains(textbox.Text.ToLower()));
                if (parent.Count() == 0)
                {
                    var list = new ObservableCollection<ModTypeForTreeView>(); // new collection containing filtered mods
                    var parents = collection.Where(p => p.Children.Any(y => y.DisplayName.ToLower().Contains(textbox.Text.ToLower())));
                    if (parents.Count() != 0)
                    {
                        foreach (var x in parents.ToList())
                        {
                            var temp = new ModTypeForTreeView(x.DisplayName, false);
                            list.Add(temp);
                            temp.Use = false;

                            // fix: use status!
                            var children = x.Children.Where(y => y.DisplayName.ToLower().Contains(textbox.Text.ToLower()));
                            foreach (var child in children)
                            {
                                temp.Children.Add(child);
                            }
                        }
                    }
                    parent = list;
                }
                tree.DataContext = parent;
            }
        }
    }
}