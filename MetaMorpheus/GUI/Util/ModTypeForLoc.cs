using System.ComponentModel;

namespace MetaMorpheusGUI
{
    public class ModTypeForLoc : INotifyPropertyChanged
    {
        private bool? _isChecked;

        public ModTypeForLoc(string displayName)
        {
            DisplayName = displayName;
        }

        public event PropertyChangedEventHandler PropertyChanged;

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

        protected void RaisePropertyChanged(string name)
        {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(name));
        }

        private void SetUseStatus(bool? value)
        {
            _isChecked = value;

            RaisePropertyChanged("Use");
        }
    }
}