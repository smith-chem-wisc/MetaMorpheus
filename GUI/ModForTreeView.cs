using System.ComponentModel;
using System.Windows.Media;

namespace MetaMorpheusGUI
{
    public class ModForTreeView : INotifyPropertyChanged
    {
        private bool _IsChecked;

        public ModForTreeView(string toolTip, bool use, string displayName, bool bad, ModTypeForTreeView parent)
        {
            ToolTipStuff = toolTip;
            Parent = parent;
            Use = use;
            DisplayName = displayName;
            if (bad)
                Background = new SolidColorBrush(Colors.Red);
            else
                Background = new SolidColorBrush(Colors.Transparent);
        }

        public event PropertyChangedEventHandler PropertyChanged;

        public ModTypeForTreeView Parent { get; }
        public string ToolTipStuff { get; }

        public bool Use
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
        public Brush Background { get; }

        internal void SetUseStatus(bool value)
        {
            if (value == Use)
                return;

            _IsChecked = value;
            Parent.VerifyCheckState();

            RaisePropertyChanged("Use");
        }

        protected void RaisePropertyChanged(string name)
        {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(name));
        }
    }
}