using System.ComponentModel;
using System.Linq;
using System.Windows.Media;

namespace MetaMorpheusGUI
{
    public class ModForTreeView : INotifyPropertyChanged
    {
        private bool _isChecked;

        public ModForTreeView(string toolTip, bool use, string modName, bool bad, ModTypeForTreeView parent)
        {
            ToolTipStuff = toolTip;
            Parent = parent;
            Use = use;
            ModName = modName;

            DisplayName = modName;

            if (toolTip.ToLower().Contains("terminal"))
            {
                var split = toolTip.Split('\n');

                string location = split.First(p => p.StartsWith("PP   "));
                location = location.Substring(5, location.Length - 5).Trim();

                switch (location)
                {
                    case "N-terminal.": DisplayName += " (Prot N-Term)"; break;
                    case "C-terminal.": DisplayName += " (Prot C-Term)"; break;
                    case "Peptide N-terminal.": DisplayName += " (Pep N-Term)"; break;
                    case "Peptide C-terminal.": DisplayName += " (Pep C-Term)"; break;
                }
            }

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
                return _isChecked;
            }
            set
            {
                SetUseStatus(value);
            }
        }

        public string ModName { get; }
        public string DisplayName { get; }
        public Brush Background { get; }

        internal void SetUseStatus(bool value)
        {
            if (value == Use)
                return;

            _isChecked = value;
            Parent.VerifyCheckState();

            RaisePropertyChanged("Use");
        }

        protected void RaisePropertyChanged(string name)
        {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(name));
        }
    }
}