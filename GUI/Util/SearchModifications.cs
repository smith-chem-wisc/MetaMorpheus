using System;
using System.Collections.ObjectModel;
using System.Linq;
using System.Windows.Controls;
using System.Windows.Threading;

namespace MetaMorpheusGUI.Util
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
                    var l = new ObservableCollection<ModForTreeView>();
                    foreach (var d in collection) // lists all child mods with user-specified keyword
                    {
                        var list = d.Children.Where(x => x.DisplayName.ToLower().Contains(textbox.Text.ToLower()));
                        if (list.Count() != 0)
                        {
                            list.ToList().ForEach(x => l.Add(x));
                        }
                    }
                    tree.DataContext = l;
                }
                else
                {
                    tree.DataContext = parent; // shows parent node (not expanded)
                }
            }
        }
    }
}