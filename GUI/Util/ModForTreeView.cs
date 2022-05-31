using GuiFunctions;
using OxyPlot;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Linq;
using System.Windows.Media;

namespace MetaMorpheusGUI
{
    public class ModForTreeView : BaseViewModel
    {
        #region Private Properties

        private bool _isChecked;
        private string _selectedColor;
        private SolidColorBrush _colorBrush;

        #endregion

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
                _isChecked = value;
                OnPropertyChanged(nameof(Use));
            }
        }

        public string ModName { get; }
        public bool HasChanged { get; set; }
        public string DisplayName { get; }
        public Brush Background { get; }
        public ObservableCollection<string> PossibleColors { get; set; }
        public string SelectedColor
        {
            get { return _selectedColor; }
            set
            {
                _selectedColor = value;
                ColorBrush = DrawnSequence.ParseColorBrushFromName(_selectedColor);
                OnPropertyChanged(nameof(SelectedColor));
            }
        }
        public SolidColorBrush ColorBrush
        {
            get { return _colorBrush; }
            set
            {
                _colorBrush = value;
                OnPropertyChanged(nameof(ColorBrush));
            }
        }

        #endregion

        #region Constructor

        /// <summary>
        /// Constructor for use in GPTMD task window
        /// </summary>
        /// <param name="toolTip"></param>
        /// <param name="use"></param>
        /// <param name="modName"></param>
        /// <param name="bad"></param>
        /// <param name="parent"></param>
        public ModForTreeView(string toolTip, bool use, string modName, bool bad, ModTypeForTreeView parent, ObservableCollection<string> colors = null)
        {
            ToolTipStuff = toolTip;
            Parent = parent;
            Use = use;
            ModName = modName;
            DisplayName = modName;
            PossibleColors = colors;
            AddSpaces(PossibleColors);
            OxyColor color = MetaDrawSettings.ModificationTypeToColor[modName];
            SelectedColor = AddSpaces(color.GetColorName());
            ColorBrush = DrawnSequence.ParseColorBrushFromOxyColor(color);

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

        #endregion

        public void SelectionChanged(string newColor)
        {
            SelectedColor = newColor;
            HasChanged = true;
        }
    }
}