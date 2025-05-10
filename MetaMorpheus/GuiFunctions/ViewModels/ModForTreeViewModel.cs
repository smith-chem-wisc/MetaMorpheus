using GuiFunctions;
using OxyPlot;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Linq;
using System.Windows.Media;

namespace GuiFunctions
{
    public class ModForTreeViewModel : BaseViewModel
    {
        #region Private Properties

        private bool _isChecked;
        private string _selectedColor;
        private SolidColorBrush _colorBrush;
        private bool _colorDataLoaded;

        #endregion

        #region Public Properties
        public ModTypeForTreeViewModel Parent { get; }
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
        public string SelectedColor
        {
            get
            {
                EnsureColorDataLoaded();
                return _selectedColor;
            }
            set
            {
                _selectedColor = value;
                ColorBrush = DrawnSequence.ParseColorBrushFromName(_selectedColor);
                _colorDataLoaded = true;
                OnPropertyChanged(nameof(SelectedColor));
            }
        }
        public SolidColorBrush ColorBrush
        {
            get
            {
                EnsureColorDataLoaded();
                return _colorBrush;
            }
            set
            {
                _colorBrush = value;
                OnPropertyChanged(nameof(ColorBrush));
            }
        }

        #endregion

        #region Constructor

        public ModForTreeViewModel(string toolTip, bool use, string modName, bool bad, ModTypeForTreeViewModel parent)
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

        #endregion

        /// <summary>
        /// Ensures the color data is loaded. This is necessary because the color data is not loaded until the first time it is accessed.
        /// This enables the use of the same control in MetaDraw and Task Windows without loading the color data for task windows. 
        /// </summary>
        private void EnsureColorDataLoaded()
        {
            if (!_colorDataLoaded)
            {
                if (MetaDrawSettings.ModificationTypeToColor != null)
                {
                    // This if statement prevents a crash from loading a search task modifications not found on launch
                    // This can occur due to new custom modifications or a mod in the xml database that was not in our initial list
                    if (!MetaDrawSettings.ModificationTypeToColor.TryGetValue(ModName, out OxyColor color))
                        color = MetaDrawSettings.FallbackColor;
                    _selectedColor = AddSpaces(MetaDrawSettings.PossibleColors[color]);
                    _colorBrush = DrawnSequence.ParseColorBrushFromOxyColor(color);
                }
                _colorDataLoaded = true;
            }
        }

        public void SelectionChanged(string newColor)
        {
            SelectedColor = newColor;
            HasChanged = true;
        }
    }
}