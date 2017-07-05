using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Windows.Media;

namespace MetaMorpheusGUI
{
    public class ModTypeForTreeView : INotifyPropertyChanged
    {
        #region Private Fields

        private bool? _isChecked;

        #endregion Private Fields

        #region Public Constructors

        public ModTypeForTreeView(string displayName, bool bad)
        {
            Children = new ObservableCollection<ModForTreeView>();
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

        public bool? Use
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

        public ObservableCollection<ModForTreeView> Children { get; }

        public Brush Background { get; }

        #endregion Public Properties

        #region Public Methods

        public void VerifyCheckState()
        {
            bool? state = null;
            for (int i = 0; i < Children.Count; ++i)
            {
                bool current = Children[i].Use;
                if (i == 0)
                    state = current;
                else if (state != current)
                {
                    state = null;
                    break;
                }
            }
            SetUseStatus(state);
        }

        #endregion Public Methods

        #region Protected Methods

        protected void RaisePropertyChanged(string name)
        {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(name));
        }

        #endregion Protected Methods

        #region Private Methods

        private void SetUseStatus(bool? value)
        {
            if (value == Use)
                return;

            _isChecked = value;

            if (value.HasValue)
                foreach (var child in Children)
                    child.SetUseStatus(value.Value);

            RaisePropertyChanged("Use");
        }

        #endregion Private Methods
    }
}