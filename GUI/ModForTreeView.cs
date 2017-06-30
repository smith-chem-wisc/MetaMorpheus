using System.ComponentModel;
using System.Windows.Media;

namespace MetaMorpheusGUI
{
    public class ModForTreeView : INotifyPropertyChanged
    {
        #region Private Fields

        private bool _isChecked;

        #endregion Private Fields

        #region Public Constructors

        public ModForTreeView(string toolTip, bool use, string displayName, bool bad, ModTypeForTreeView _parent)
        {
            ToolTipStuff = toolTip;
            Parent = _parent;
            Use = use;
            DisplayName = displayName;
            if (bad)
                Background = new SolidColorBrush(Colors.Red);
            else
                Background = new SolidColorBrush(Colors.Transparent);
        }

        #endregion Public Constructors

        #region Public Events

        public event PropertyChangedEventHandler PropertyChanged;

        #endregion Public Events

        #region Public Properties

        public ModTypeForTreeView Parent { get; }
        public string ToolTipStuff { get; }

        public bool Use
        {
            get
            {
                return _isChecked;
            }
            set
            {
                SetUseStatus(value);
            }
        }

        public string DisplayName { get; }
        public Brush Background { get; }

        #endregion Public Properties

        #region Internal Methods

        internal void SetUseStatus(bool value)
        {
            if (value == Use)
                return;

            _isChecked = value;
            Parent.VerifyCheckState();

            RaisePropertyChanged("Use");
        }

        #endregion Internal Methods

        #region Protected Methods

        protected void RaisePropertyChanged(string name)
        {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(name));
        }

        #endregion Protected Methods
    }
}