using System.ComponentModel;

namespace MetaMorpheusGUI
{
    public class ModTypeForLoc : INotifyPropertyChanged
    {
        #region Private Fields

        private bool? _isChecked;

        #endregion Private Fields

        #region Public Constructors

        public ModTypeForLoc(string displayName)
        {
            DisplayName = displayName;
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

        #endregion Public Properties

        #region Protected Methods

        protected void RaisePropertyChanged(string name)
        {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(name));
        }

        #endregion Protected Methods

        #region Private Methods

        private void SetUseStatus(bool? value)
        {
            _isChecked = value;

            RaisePropertyChanged("Use");
        }

        #endregion Private Methods
    }
}