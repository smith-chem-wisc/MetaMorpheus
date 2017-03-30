using System.Windows.Media;

namespace MetaMorpheusGUI
{
    public class ModForTreeView
    {

        #region Public Constructors

        public ModForTreeView(bool use, string displayName, bool bad)
        {
            Use = use;
            DisplayName = displayName;
            if (bad)
                Background = new SolidColorBrush(Colors.Red);
            else
                Background = new SolidColorBrush(Colors.Transparent);
        }

        #endregion Public Constructors

        #region Public Properties

        public bool Use { get; set; }
        public string DisplayName { get; }
        public Brush Background { get; }

        #endregion Public Properties

    }
}