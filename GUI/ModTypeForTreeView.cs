using System.Collections.ObjectModel;
using System.Windows.Media;

namespace MetaMorpheusGUI
{
    public class ModTypeForTreeView
    {

        #region Public Constructors

        public ModTypeForTreeView(bool use, string displayName, bool bad)
        {
            Use = use;
            DisplayName = displayName;
            Children = new ObservableCollection<ModForTreeView>();
            if (bad)
                Background = new SolidColorBrush(Colors.Red);
            else
                Background = new SolidColorBrush(Colors.Transparent);
        }

        #endregion Public Constructors

        #region Public Properties

        public bool? Use { get; set; }
        public string DisplayName { get; }
        public ObservableCollection<ModForTreeView> Children { get; }
        public Brush Background { get; }

        #endregion Public Properties

    }
}