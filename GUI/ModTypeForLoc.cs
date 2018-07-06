using System.ComponentModel;

namespace MetaMorpheusGUI
{
    public class ModTypeForLoc : INotifyPropertyChanged
    {
        private bool? _IsChecked;

        public ModTypeForLoc(string displayName)
        {
            DisplayName = displayName;
        }

        public event PropertyChangedEventHandler PropertyChanged;

        public bool? Use
        {
            get
            {
                return _IsChecked;
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
            _IsChecked = value;
            RaisePropertyChanged("Use");
        }
    }
}