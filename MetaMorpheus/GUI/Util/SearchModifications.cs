using System;
using System.Collections.ObjectModel;
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
            FilterTree(textbox, tree, collection, null);
        }

        /// <summary>
        /// As above, but matching text of the caller's choosing instead of only
        /// <see cref="ModForTreeViewModel.ModName"/> -- see <see cref="ModTreeFilter.Filter"/>.
        /// </summary>
        public static void FilterTree(TextBox textbox, TreeView tree, ObservableCollection<ModTypeForTreeViewModel> collection,
            Func<ModForTreeViewModel, string> searchText)
        {
            string key = textbox.Text.ToLower();
            if (string.IsNullOrEmpty(key))
            {
                tree.DataContext = collection; // shows full tree if nothing is searched
                return;
            }

            tree.DataContext = ModTreeFilter.Filter(collection, key, searchText);
        }
    }
}
