using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Windows.Media;

namespace MetaMorpheusGUI
{
    public class ModTypeForTreeView : INotifyPropertyChanged
    {
        private bool? _isChecked;

        public ModTypeForTreeView(string displayName, bool bad)
        {
            Children = new ObservableCollection<ModForTreeView>();
            Expanded = false;
            DisplayName = displayName;
            if (bad)
                Background = new SolidColorBrush(Colors.Red);
            else
                Background = new SolidColorBrush(Colors.Transparent);
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

        public bool Expanded { get; set; }

        public string DisplayName { get; }

        public ObservableCollection<ModForTreeView> Children { get; }

        public Brush Background { get; }

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

        protected void RaisePropertyChanged(string name)
        {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(name));
        }

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
    }
}